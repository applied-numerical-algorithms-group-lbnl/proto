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
    m_stensizes.resize(a_dstVoFs.size());

    Box srcbox = a_validBox.grow(a_srcGhost);
    Box dstbox = a_validBox.grow(a_dstGhost);

    EBBoxData<srcCenter, data_t, 1> srcData(srcbox, a_graph);
    EBBoxData<srcCenter, data_t, 1> dstData(dstbox, a_graph);
#ifdef PROTO_CUDA
    vector<pair_t> h_serializedStenInfo;
#endif
    for (int idst = 0; idst < a_dstVoFs.size(); idst++)
    {
      m_stensizes[idst] = a_stencil[idst].size();
      m_dstAccess[idst] = dstData.indexInfo(a_dstVoFs[idst], 0);
      m_ebstencil[idst].resize(a_stencil[idst].size());

      const LocalStencil<srcCenter, data_t>& stenLoc = a_stencil[idst];
      for (int isten = 0; isten < stenLoc.size(); isten++)
      {
        auto entry = stenLoc.m_entries[isten];
        m_ebstencil[idst][isten].first  = srcData.indexInfo(entry.m_vof, 0);
        m_ebstencil[idst][isten].second = entry.m_weight;
#ifdef PROTO_CUDA
        h_serializedStenInfo.push_back(m_ebstencil[idst][isten]);
#endif
      }
    }
#ifdef PROTO_CUDA
    cudaPrep(h_serializedStenInfo);
#endif
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

#ifdef PROTO_CUDA
    cudaApply(a_lph, a_phi, a_scale);
#else
    hostApply(a_lph, a_phi, a_scale);
#endif
  }

///
#ifdef PROTO_CUDA
  template <CENTERING srcCenter, CENTERING dstCenter, typename data_t>
  void 
  AggStencil<srcCenter, dstCenter, data_t>::
  cudaPrep(const vector<pair_t>& h_serializedStenInfo)

  {
    uint64_t a_numStenEntriesTotal = h_serializedStenInfo.size();
    
    uint64_t totStenSize = a_numStenEntriesTotal*(sizeof(pair_t));
    uint64_t destAccSize = numptsDst*(sizeof(EBDataLoc));
    uint64_t stenLenSize = numptsDst*(sizeof(uint64_t));
    ///allocate memory on the device
    cudaMalloc(&d_ebstencil, totStenSize);
    cudaMalloc(&d_dstAccess, destAccSize);
    cudaMalloc(&d_stensizes, stenLenSize);

    ///copy from the host
    cudaMemcpy(d_dstAccess, m_dstAccess.data(),           destAccSize, cudaMemcpyHostToDevice);
    cudaMemcpy(d_stensizes, m_stensizes.data(),           stenLenSize, cudaMemcpyHostToDevice);
    cudaMemcpy(d_ebstencil, h_serializedStenInfo.data(),  totStenSize, cudaMemcpyHostToDevice);

    //figure out blocks and strides
    m_stride = 32;
    bool found = false;
    while(!found)
    {
      if((numPtsDst % m_stride) == 0)
      {
        found = true;
      }
      else
      {
        m_stride /= 2;
      }
    }
    m_blocks = numPtsDst/m_stride;
  }
/////

  template <CENTERING srcCenter, CENTERING dstCenter, typename data_t>
  void 
  AggStencil<srcCenter, dstCenter, data_t>::
  cudaRelease()
  {
    cudaFree(d_ebstencil);
    cudaFree(d_dstAccess);
    cudaFree(d_stensizes);
  }
/////

  template <CENTERING srcCenter, CENTERING dstCenter, typename data_t>
  void 
  AggStencil<srcCenter, dstCenter, data_t>::
  cudaApply(EBBoxData<      dstCenter, data_t, 1>     & a_lph,
            const EBBoxData<srcCenter, data_t, 1>     & a_phi,
            const data_t                              & a_scale)
  {
  }
#endif
}
///////////
#endif
