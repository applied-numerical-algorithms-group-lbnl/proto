
#ifndef _Proto_EBStencilImplem_H_
#define _Proto_EBStencilImplem_H_
#include "Proto_EBHostData.H"
namespace Proto
{
  ///
  /**
     regular stencil gets applied in the interior of the domain.
     The aggstencil overwrites it near the EB.
     The domainBC stencils fill ghost cells outside the domain.  
     The vector length needs to be the same as the one-sided span of the stencil.
  */
  template <int order, typename data_t, CENTERING srcCenter, CENTERING dstCenter>
  EBStencil<order, data_t, srcCenter, dstCenter>::
  EBStencil(string              a_stencilName,
            string              a_domainBCName,
            string              a_ebbcName,
            const Box         & a_validBox,
            const Box         & a_domain,
            const Point       & a_srcGhost,
            const Point       & a_dstGhost,
            const EBGraph     & a_graph,
            const VoluData    & a_voludata,
            const EBFaData    & a_ebfadata,
            const XFacData    & a_xfacdata,
            const YFacData    & a_yfacdata,
            const ZFacData    & a_zfacdata,
            double              a_dx,
            bool                a_periodiconly,
            bool                a_needDiagonalWeights)
  {
    m_periodic_domain_bcs_only = a_periodiconly;
    vector< EBIndex<dstCenter> >              dstVoFs;
    vector< LocalStencil<srcCenter, data_t> > stencil;

    m_hasDomainBCs = 
      archive_t::getStencil(dstVoFs, stencil, m_regStencilInterior, 
                            m_regStencilBCS, m_BCApplyBoxes, m_regApplyBox,
                            a_stencilName,  a_domainBCName, a_ebbcName,    
                            a_validBox, a_domain,      
                            a_srcGhost, a_dstGhost,    
                            a_graph, a_voludata, a_ebfadata,    
                            a_xfacdata, a_yfacdata, a_zfacdata,    
                            a_dx, a_periodiconly);

    m_irregStencil = shared_ptr<AggStencil<srcCenter, dstCenter, data_t> >
      (new AggStencil<srcCenter, dstCenter, data_t>(dstVoFs, stencil, a_graph,
                                                    a_validBox, a_srcGhost, a_dstGhost));

    if(a_needDiagonalWeights)
    {
      computeDiagonalWeights(dstVoFs, stencil, a_graph, a_validBox);
    }
  }
  ///
  template <int order, typename data_t, CENTERING srcCenter, CENTERING dstCenter>
  void 
  EBStencil<order, data_t, srcCenter, dstCenter>::
  computeDiagonalWeights(const vector< EBIndex<dstCenter> >              & a_dstVoFs,
                         const vector< LocalStencil<srcCenter, data_t> > & a_stencil,
                         const EBGraph                                   & a_graph,
                         const Box                                       & a_validBox)
  {
    //not sure this makes any sense if the centerings are different
    if(srcCenter == dstCenter)
    {
      m_diagonalWeights = shared_ptr< EBBoxData<dstCenter, data_t, 1> >
        (new EBBoxData<dstCenter, data_t, 1>(a_validBox, a_graph));
    
      //have to compute stuff on the host first
      EBHostData<dstCenter, data_t, 1> hostDiagonal(a_validBox, a_graph);
      //first set the regular data to the value of the regular stencil's diagonal.
      data_t regDiagonal = m_regStencilInterior.diagonalValue();
      hostDiagonal.setVal(regDiagonal);
    
      //now fix up values on irregular cells
      for(int idst = 0; idst < a_dstVoFs.size(); idst++)
      {
        const EBIndex<dstCenter>& dstvof               = a_dstVoFs[idst];
        const LocalStencil<srcCenter, data_t>& irrsten = a_stencil[idst];
        data_t diagval = 0;
        for(int ientry = 0; ientry < irrsten.size(); ientry++)
        {
          const EBIndex<srcCenter> srcind = irrsten.m_entries[ientry].m_vof;
          EBIndex<dstCenter> dststen = getEBIndex<dstCenter, srcCenter>(srcind);
          if(dststen == dstvof)
          {
            diagval += irrsten.m_entries[ientry].m_weight;
          }
        }
        hostDiagonal(dstvof, 0) = diagval;
      }

      //now copy the diagonal value to the device
      EBHostData<dstCenter, data_t, 1>& hostfab = hostDiagonal;
      EBBoxData< dstCenter, data_t, 1>& devifab = *m_diagonalWeights;
      for(unsigned int iwhichfab = 0; iwhichfab < 2; iwhichfab++)
      {
        size_t   fabsizeHost, fabsizeDevi; 
        data_t*  dataptrHost, *dataptrDevi;
        hostfab.dataInfo(dataptrHost, fabsizeHost, iwhichfab);
        devifab.dataInfo(dataptrDevi, fabsizeDevi, iwhichfab);

        size_t datSize = sizeof(data_t)*fabsizeDevi;

#ifdef PROTO_CUDA
        cudaMemcpy(dataptrDevi, dataptrHost, datSize, cudaMemcpyHostToDevice);
        cudaError err = cudaGetLastError();
        if (err != cudaSuccess)
        {
          fprintf(stderr, "cudaGetLastError() failed at %s:%i : %s\n",
                  __FILE__, __LINE__, cudaGetErrorString(err));
        }
#else
        memcpy(dataptrDevi, dataptrHost, datSize);
#endif  
      }
    }
  }

  ///
  template <int order, typename data_t, CENTERING srcCenter, CENTERING dstCenter>
  void 
  EBStencil<order, data_t, srcCenter, dstCenter>::
  apply(EBBoxData<       dstCenter, data_t, 1>  &   a_dst,
        const EBBoxData< srcCenter, data_t, 1>  &   a_src,
        bool            a_initToZero,
        double          a_scale)  const
  {
    PR_TIME("EBStencil::apply");
    if(a_initToZero)
    {
      a_dst.setVal(0.);
    }
    
    //first set the boundary conditions at the domain box
    //this fills the ghost cells when the box abuts the domain.
    if(m_hasDomainBCs)
    {
      BoxData<data_t, 1>& castSrc = const_cast<BoxData<data_t, 1> &>(a_src.getRegData());
      for(int ibc = 0; ibc < 2*DIM; ibc++)
      {
        for(int ivec = 0; ivec < m_regStencilBCS[ibc].size(); ivec++)
        {
          Box bcsrcbox = m_BCApplyBoxes[ibc][ivec];
          m_regStencilBCS[ibc][ivec].apply(castSrc, castSrc,
                                           bcsrcbox, true, 1.0);
        }
      }
    }
    bool incrementOnly = !a_initToZero;
    //already set stuff to zero if necessary
    m_regStencilInterior.apply(a_src.getRegData(), 
                               a_dst.getRegData(),
                               m_regApplyBox, incrementOnly, a_scale);

    m_irregStencil->apply(a_dst, a_src, a_scale, incrementOnly);
    
  }
              
}

#endif
