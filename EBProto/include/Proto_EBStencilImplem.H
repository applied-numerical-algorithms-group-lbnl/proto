
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
  EBStencil<typename data_t, CENTERING srcCenter, CENTERING dstCenter>::
  EBStencil(string              a_stencilName,
            string              a_domainBCName[2*DIM],
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
            bool periodiconly)
  {
    HERE!!!
  }

  ///
  void 
  EBStencil<typename data_t, CENTERING srcCenter, CENTERING dstCenter>::
  apply(EBBoxData<       dstCenter, data_t, 1> > &   a_dst,
        const EBBoxData< srcCenter, data_t, 1> > &   a_src,
        bool            a_initToZero,
        const double    a_scale)  const
  {
  }
               
              
}

#endif
