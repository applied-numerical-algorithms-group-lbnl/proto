#ifndef __Proto_AggStencilImplem__
#define __Proto_AggStencilImplem__

namespace Proto
{

  ///
  template <CENTERING srcCenter, CENTERING dstCenter, typename data_t>
  AggStencil<srcCenter, dstCenter, data_t>::
  AggStencil(const vector< EBIndex<     dstCenter        > > & a_dstVoFs,
             const vector< LocalStencil<srcCenter, data_t> > & a_stencil,
             const EBGraph                                   & a_graph,
             const Box                                       & a_validBox,
             const Point                                     & a_srcGhost,
             const Point                                     & a_dstGhost)
  {
    PR_TIME("AggSten.constructor");
    m_ebstencil.resize(a_dstVoFs.size());
    m_dstAccess.resize(a_dstVoFs.size());

    Box srcbox = a_validBox.grow(a_srcGhost);
    Box dstbox = a_validBox.grow(a_dstGhost);

    EBBoxData<srcCenter, data_t, 1> srcData(srcbox, a_graph);
    EBBoxData<srcCenter, data_t, 1> dstData(dstbox, a_graph);

    for (int idst = 0; idst < a_dstVoFs.size(); idst++)
    {
      m_dstAccess[idst] = dstData.indexInfo(a_dstVoFs[idst], 0);
      m_ebstencil[idst].resize(a_stencil[idst].size());

      const LocalStencil<srcCenter, data_t>& stenLoc = a_stencil[idst];
      for (int isten = 0; isten < stenLoc.size(); isten++)
      {
        auto entry = stenLoc.m_entries[isten];
        m_ebstencil[idst][isten].first  = srcData.indexInfo(entry.m_vof, 0);
        m_ebstencil[idst][isten].second = entry.m_weight;
      }
    }
  }


  /**
     lphi = L(phi)
  */
  template <CENTERING srcCenter, CENTERING dstCenter, typename data_t>
  void 
  AggStencil<srcCenter, dstCenter, data_t>::
  apply(EBBoxData<      dstCenter, data_t, 1>     & a_lph,
        const EBBoxData<srcCenter, data_t, 1>     & a_phi,
        const data_t                              & a_scale)
  {
    PR_TIME("AggSten::apply");

    data_t*       dataPtrsLph[2];
    const data_t* dataPtrsPhi[2];

    int varDst = 0; int varSrc = 0; //we might have to revisit this
    for (int ivec = 0; ivec < 2; ivec++)
    {
      dataPtrsLph[ivec] = a_lph.dataPtr(ivec, varDst);
      dataPtrsPhi[ivec] = a_phi.dataPtr(ivec, varSrc);
    }

    for (int idst = 0; idst < m_ebstencil.size(); idst++)
    {
      data_t* lphiPtr =  dataPtrsLph[m_dstAccess[idst].m_dataID] + m_dstAccess[idst].m_offset;
      data_t& lphi = *lphiPtr;

      //initializing to zero is necessary here.   incrementOnly is not an option right now
      lphi =  0.;

      const stencil_t& ebstencil = m_ebstencil[idst];
      for (int isten = 0; isten < ebstencil.size(); isten++)
      {
        const data_t & weight  = ebstencil[isten].second;
        const size_t & offset  = ebstencil[isten].first.m_offset;
        const int    & dataID  = ebstencil[isten].first.m_dataID;
        const data_t  & phiVal = *(dataPtrsPhi[dataID] + offset);
        lphi += phiVal*weight;
      }
      PR_FLOPS(ebstencil.size()*2);
    }
  }
}
///////////
#endif
