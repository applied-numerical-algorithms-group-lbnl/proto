
#ifndef _Proto_EBStencilImplem_H_
#define _Proto_EBStencilImplem_H_

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
            bool                a_periodiconly)
  {
    m_periodic_domain_bcs_only = a_periodiconly;
    vector<EBIndex<dstCenter> > dstVoFs;
    vector<LocalStencil<srcCenter, data_t> > stencil;

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
    
    BoxData<data_t, 1>& castSrc = const_cast<BoxData<data_t, 1> &>(a_src.getRegData());
    //first set the boundary conditions at the domain box
    //this fills the ghost cells when the box abuts the domain.
    for(int ibc = 0; ibc < 2*DIM; ibc++)
    {
      for(int ivec = 0; ivec < m_regStencilBCS[ibc].size(); ivec++)
      {
        Box bcsrcbox = m_BCApplyBoxes[ibc][ivec];
        m_regStencilBCS[ibc][ivec].apply(castSrc, castSrc,
                                 bcsrcbox, false, 1.0);
      }
    }

    bool incrementOnly = !a_initToZero;
    //already set stuff to zero if necessary
    m_regStencilInterior.apply(a_src.getRegData(), 
                               a_dst.getRegData(),
                               m_regApplyBox, incrementOnly, a_scale);

    m_irregStencil->apply(a_dst, a_src, incrementOnly, a_scale);
    
  }
               
              
}

#endif
