
#ifndef _Proto_EBStencilArchiveImplem_H_
#define _Proto_EBStencilArchiveImplem_H_

namespace Proto
{
  inline void 
  get2dOrderPoissonStencil(vector<EBIndex<CELL> >     & a_dstVoFs,
                           vector<LocalStencil<CELL> >& a_stencil,
                           Stencil<double>            & a_regStencilInterior,
                           Box                        & a_regApplyBox,
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
    if(a_periodiconly)
    {
      a_regApplyBox = a_validBox;
    }
    else
    {
      PR_error("domain bcs not implemented here yet");
    }
    a_regStencilInterior =   (1.0/(a_dx*a_dx))*(Stencil<double>::Laplacian());

    a_dstVoFs = a_graph.getIrregVoFs(a_validBox);
    a_stencil.resize(a_dstVoFs.size());
    for(int ivof = 0; ivof < a_dstVoFs.size(); ivof++)
    {
      HERE!!
    }
  }
  ///this deep down in the bowels, everything is a template specialization
  template < >
  void
  EBStencilArchive<CELL, CELL>::
  getStencil(vector<EBIndex<dstCenter> >     & a_dstVoFs,
             vector<LocalStencil<srcCenter> >& a_stencil,
             Stencil<double>                 & a_regStencilInterior,
             vector<Stencil<double> >          a_regStencilBCS[2*DIM],
             Box                             & a_regApplyBox,
             string              a_stencilName,
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
             bool                a_periodiconly)
  {
    if(a_stencilName == string("Second_Order_Poisson"))
    {
      get2dOrderPoissonStencil(a_dstVoFs, a_stencil, a_regStencilInterior, a_applyBox,
                               a_ebbcName,
                               a_validBox,
                               a_domain,
                               a_srcGhost,
                               a_dstGhost,
                               a_graph,
                               a_voludata,
                               a_ebfadata,
                               a_xfacdata,
                               a_yfacdata,
                               a_zfacdata,
                               a_periodiconly);
    }
    else
    {
      PR_error("stencil_archive_cell_cell::stencil string not found");
    }
  }


}

#endif
