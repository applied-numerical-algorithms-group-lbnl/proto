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

    //begin debugging hooks
    m_dstVoFs = a_dstVoFs;
    m_stencil = a_stencil;
    //end   debugging hooks

    Box srcbox = a_validBox.grow(a_srcGhost);
    Box dstbox = a_validBox.grow(a_dstGhost);

    m_flopsThisApply = 0;        

    EBBoxData<srcCenter, data_t, 1> srcData(srcbox, a_graph);
    EBBoxData<srcCenter, data_t, 1> dstData(dstbox, a_graph);

    vector<uint64_t>  h_stenSizes(a_dstVoFs.size());
    vector<uint64_t>  h_stenStart(a_dstVoFs.size());
    vector<pair_t<data_t> >    h_serializedStenInfo;
    uint64_t startLoc = 0;

    for (int idst = 0; idst < a_dstVoFs.size(); idst++)
    {
      m_dstAccess[idst] = dstData.indexInfo(a_dstVoFs[idst], 0);
      m_ebstencil[idst].resize(a_stencil[idst].size());

      const LocalStencil<srcCenter, data_t>& stenLoc = a_stencil[idst];
      m_flopsThisApply += 3*stenLoc.size();
//      m_flopsThisApply += 10;
      data_t sumweight = 0;
      for (int isten = 0; isten < stenLoc.size(); isten++)
      {
        auto entry = stenLoc.m_entries[isten];
        m_ebstencil[idst][isten].first  = srcData.indexInfo(entry.m_vof, 0);
        m_ebstencil[idst][isten].second = entry.m_weight;
        sumweight += entry.m_weight;
        h_serializedStenInfo.push_back(m_ebstencil[idst][isten]);
      }
//      printf("agd:sumweight = %f\n", sumweight);
      

      h_stenSizes[idst] = a_stencil[idst].size();
      h_stenStart[idst] = startLoc;
      startLoc         += a_stencil[idst].size();

    }
#ifdef PROTO_CUDA
    cudaPrep(h_serializedStenInfo, h_stenSizes, h_stenStart);
#else
    hostPrep(h_serializedStenInfo, h_stenSizes, h_stenStart);
#endif
//    printf("dstvofs size = %d, flops per apply = %d\n ",(int)a_dstVoFs.size(), m_flopsThisApply);
  }


  template <CENTERING srcCenter, CENTERING dstCenter, typename data_t>
  void 
  AggStencil<srcCenter, dstCenter, data_t>::
  hostPrep(const vector<pair_t<data_t> >   &   h_serializedStenInfo,
           const vector<uint64_t>          &   h_stenSizes,
           const vector<uint64_t>          &   h_stenStart)

  {
    uint64_t a_numStenEntriesTotal = h_serializedStenInfo.size();
    uint64_t numptsDst = m_dstAccess.size();
    uint64_t totStenSize = a_numStenEntriesTotal*(sizeof(pair_t<data_t>));
    uint64_t destAccSize = numptsDst*(sizeof(EBDataLoc));
    uint64_t stenLenSize = numptsDst*(sizeof(uint64_t));
    uint64_t stStartSize = numptsDst*(sizeof(uint64_t));
    ///allocate memory on the device
    d_ebstencil = (pair_t<data_t>*)malloc(totStenSize);
    d_dstaccess = (EBDataLoc*     )malloc(destAccSize);
    d_stensizes = (uint64_t*      )malloc(stenLenSize);
    d_stenstart = (uint64_t*      )malloc(stStartSize);

    ///copy from the host
    memcpy(d_dstaccess, m_dstAccess.data(),           destAccSize);
    memcpy(d_stensizes, h_stenSizes.data(),           stenLenSize);
    memcpy(d_stenstart, h_stenStart.data(),           stStartSize);
    memcpy(d_ebstencil, h_serializedStenInfo.data(),  totStenSize);

  }
  /**
     lphi = L(phi)
  */
  template <CENTERING srcCenter, CENTERING dstCenter, typename data_t>
  void 
  AggStencil<srcCenter, dstCenter, data_t>::
  apply(EBBoxData<      dstCenter, data_t, 1>     & a_lph,
        const EBBoxData<srcCenter, data_t, 1>     & a_phi,
        const data_t                              & a_scale,
        bool                                        a_incrementOnly) const
  {
    PR_TIME("AggSten::apply");

#ifdef PROTO_CUDA
    cudaApply(a_lph, a_phi, a_scale, a_incrementOnly);
#else
    hostApply(a_lph, a_phi, a_scale, a_incrementOnly);
#endif
    PR_FLOPS(m_flopsThisApply);
  }
///
  template <typename data_t>
  void 
  hostAggStencilIndexer(int a_begin, int a_end, int idx,
                        const data_t*           a_dataPtrsSrc[2],
                        data_t*                 a_dataPtrsDst[2],
                        const pair_t<data_t>*   a_ebstencil,
                        const uint64_t*         a_stensizes,
                        const uint64_t*         a_stenstart,
                        const EBDataLoc*        a_dstaccess,
                        bool                    a_incrementOnly)
  {
    if((idx >= a_begin) && (idx< a_end))
    {
      

      uint64_t startInd = a_stenstart[idx];
      uint64_t stenSize = a_stensizes[idx];
      const EBDataLoc& dstAcc = a_dstaccess[idx];
      data_t& lphi = a_dataPtrsDst[dstAcc.m_dataID][dstAcc.m_offset];

      if(!a_incrementOnly)
      {
        lphi = 0;
      }
      data_t sumweight = 0;
      for(uint64_t isten = 0; isten < stenSize; isten++)
      {
        const pair_t<data_t>& entry = a_ebstencil[startInd + isten];
        const EBDataLoc& srcAcc = entry.first;
        const data_t   & weight = entry.second;
        sumweight += weight;
        const data_t& phi = a_dataPtrsSrc[srcAcc.m_dataID][srcAcc.m_offset];
        lphi += weight*phi;
      }
      //printf("agi:sumweight = %f\n", sumweight);
    }
  }
  template <CENTERING srcCenter, CENTERING dstCenter, typename data_t>
  void 
  AggStencil<    srcCenter, dstCenter, data_t>::
  hostApply(EBBoxData<      dstCenter, data_t, 1>     & a_lph,
            const EBBoxData<srcCenter, data_t, 1>     & a_phi,
            const data_t                              & a_scale,
            bool                                        a_incrementOnly) const

  {
    int begin = 0;
    int end   = m_dstAccess.size();
    const data_t* dataPtrsSrc[2];
    data_t*       dataPtrsDst[2];
    for(int idat = 0; idat < 2; idat++)
    {
      dataPtrsSrc[idat]  = (a_phi.dataPtr(idat, 0));
      dataPtrsDst[idat]  = (a_lph.dataPtr(idat, 0));
    }
    for(int idx = 0; idx < end; idx++)
    {
      //begin debugging hooks
      const auto& localsten = m_stencil[idx];
      const auto& vof       = m_dstVoFs[idx];
      Point debpt(26,  17);
      int ideb = 0;
      if(vof.m_pt == debpt)
      {
        ideb = 1;
      }
      //end debugging hooks
      hostAggStencilIndexer(begin, end, idx, dataPtrsSrc, dataPtrsDst, d_ebstencil, d_stensizes, d_stenstart, d_dstaccess, a_incrementOnly);
    }
  }

  template <CENTERING srcCenter, CENTERING dstCenter, typename data_t>
  void 
  AggStencil<srcCenter, dstCenter, data_t>::
  hostRelease()
  {
    free(d_ebstencil);
    free(d_dstaccess);
    free(d_stensizes);
    free(d_stenstart);
  }
///
#ifdef PROTO_CUDA
  template <CENTERING srcCenter, CENTERING dstCenter, typename data_t>
  void 
  AggStencil<srcCenter, dstCenter, data_t>::
  cudaPrep(const vector<pair_t<data_t> >   &   h_serializedStenInfo,
           const vector<uint64_t>          &   h_stenSizes,
           const vector<uint64_t>          &   h_stenStart)

  {
    uint64_t a_numStenEntriesTotal = h_serializedStenInfo.size();
    uint64_t numptsDst = m_dstAccess.size();
    uint64_t totStenSize = a_numStenEntriesTotal*(sizeof(pair_t<data_t>));
    uint64_t destAccSize = numptsDst*(sizeof(EBDataLoc));
    uint64_t stenLenSize = numptsDst*(sizeof(uint64_t));
    uint64_t stStartSize = numptsDst*(sizeof(uint64_t));
    ///allocate memory on the device
    cudaMalloc(&d_ebstencil, totStenSize);
    cudaMalloc(&d_dstaccess, destAccSize);
    cudaMalloc(&d_stensizes, stenLenSize);
    cudaMalloc(&d_stenstart, stStartSize);

    ///copy from the host
    cudaMemcpy(d_dstaccess, m_dstAccess.data(),           destAccSize, cudaMemcpyHostToDevice);
    cudaMemcpy(d_stensizes, h_stenSizes.data(),           stenLenSize, cudaMemcpyHostToDevice);
    cudaMemcpy(d_stenstart, h_stenStart.data(),           stStartSize, cudaMemcpyHostToDevice);
    cudaMemcpy(d_ebstencil, h_serializedStenInfo.data(),  totStenSize, cudaMemcpyHostToDevice);

    //figure out blocks and strides
    m_stride = 32;
    bool found = false;
    while(!found)
    {
      if((numptsDst % m_stride) == 0)
      {
        found = true;
      }
      else
      {
        m_stride /= 2;
      }
    }
    m_blocks = numptsDst/m_stride;
  }
/////

  template <CENTERING srcCenter, CENTERING dstCenter, typename data_t>
  void 
  AggStencil<srcCenter, dstCenter, data_t>::
  cudaRelease()
  {
    cudaFree(d_ebstencil);
    cudaFree(d_dstaccess);
    cudaFree(d_stensizes);
    cudaFree(d_stenstart);
  }
/////
  template <typename data_t>
  __global__
  void 
  aggStencilIndexer(int a_begin, int a_end,
                    const data_t*           a_dataPtrRegSrc,
                    const data_t*           a_dataPtrIrrSrc,
                    data_t*                 a_dataPtrRegDst,
                    data_t*                 a_dataPtrIrrDst,
                    const pair_t<data_t>*   a_ebstencil,
                    const uint64_t*         a_stensizes,
                    const uint64_t*         a_stenstart,
                    const EBDataLoc*        a_dstaccess,
                    bool                    a_incrementOnly)
  {
    int idx = a_begin + threadIdx.x+blockIdx.x*blockDim.x;
    if((idx >= a_begin) && (idx< a_end))
    {
      const data_t* dataPtrsSrc[2];
      data_t*       dataPtrsDst[2];
      dataPtrsSrc[0]  = a_dataPtrRegSrc;
      dataPtrsDst[0]  = a_dataPtrRegDst;
      dataPtrsSrc[1]  = a_dataPtrIrrSrc;
      dataPtrsDst[1]  = a_dataPtrIrrDst;

      uint64_t startInd = a_stenstart[idx];
      uint64_t stenSize = a_stensizes[idx];
      const EBDataLoc& dstAcc = a_dstaccess[idx];

      data_t& lphi = dataPtrsDst[dstAcc.m_dataID][dstAcc.m_offset];
      if(!a_incrementOnly)
      {
        lphi = 0;
      }
      for(uint64_t isten = 0; isten < stenSize; isten++)
      {
        const pair_t<data_t>& entry = a_ebstencil[startInd + isten];
        const EBDataLoc& srcAcc = entry.first;
        const data_t   & weight = entry.second;
        const data_t& phi = dataPtrsSrc[srcAcc.m_dataID][srcAcc.m_offset];
        lphi += weight*phi;
      }
    }
  }
/////
  template <CENTERING srcCenter, CENTERING dstCenter, typename data_t>
  void 
  AggStencil<    srcCenter, dstCenter, data_t>::
  cudaApply(EBBoxData<      dstCenter, data_t, 1>     & a_lph,
            const EBBoxData<srcCenter, data_t, 1>     & a_phi,
            const data_t                              & a_scale,
            bool                    a_incrementOnly) const

  {
    int begin = 0;
    int end   = m_dstAccess.size();
    //std::cout << "cuda apply phi box = " << a_phi.getRegData().box() << std::endl;
    //std::cout << "cuda apply lph box = " << a_lph.getRegData().box() << std::endl;

    const data_t* dataPtrRegSrc = a_phi.dataPtr( 0, 0);
    const data_t* dataPtrIrrSrc = a_phi.dataPtr( 1, 0);

    data_t*       dataPtrRegDst = a_lph.dataPtr( 0, 0);
    data_t*       dataPtrIrrDst = a_lph.dataPtr( 1, 0);
    aggStencilIndexer<<<m_stride, m_blocks>>>
      (begin, end, dataPtrRegSrc, dataPtrIrrSrc, dataPtrRegDst, dataPtrIrrDst, d_ebstencil, d_stensizes, d_stenstart, d_dstaccess, a_incrementOnly);
  }
#endif
}
///////////
#endif
