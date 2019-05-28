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
#ifdef PROTO_CUDA
    vector<uint64_t>  h_stenSizes;
    vector<uint64_t>  h_stenStart;
    vector<pair_t<data_t> >    h_serializedStenInfo;
    uint64_t startLoc = 0;
#endif
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
#ifdef PROTO_CUDA
        h_serializedStenInfo.push_back(m_ebstencil[idst][isten]);
#endif
      }

#ifdef PROTO_CUDA
      h_stenSizes[idst] = a_stencil[idst].size();
      h_stenStart[idst] = a_stencil[idst].size();
      startLoc         += a_stencil[idst].size();
#endif
    }
#ifdef PROTO_CUDA
    cudaPrep(h_serializedStenInfo, h_stenSizes, h_stenStart);
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

  template <CENTERING srcCenter, CENTERING dstCenter, typename data_t>
  void 
  AggStencil<srcCenter, dstCenter, data_t>::
  hostApply(EBBoxData<      dstCenter, data_t, 1>     & a_lph,
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


///
#ifdef PROTO_CUDA
  template <CENTERING srcCenter, CENTERING dstCenter, typename data_t>
  void 
  AggStencil<srcCenter, dstCenter, data_t>::
  cudaPrep(const vector<pair_t<data_t> >   &   h_serializedStenInfo,
           const vector<uint64_t> &   h_stenSizes,
           const vector<uint64_t> &   h_stenStart)

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
                    const data_t*           a_dataPtrsSrc[2],
                    data_t*                 a_dataPtrsDst[2],
                    const pair_t<data_t>*   a_ebstencil,
                    const uint64_t*         a_stensizes,
                    const uint64_t*         a_stenstart,
                    const EBDataLoc*        a_dstaccess)
  {
    int idx = a_begin + threadIdx.x+blockIdx.x*blockDim.x;
    if((idx >= a_begin) && (idx< a_end))
    {
      uint64_t startInd = a_stenstart[idx];
      uint64_t stenSize = a_stensizes[idx];
      const EBDataLoc& dstAcc = a_dstaccess[idx];
      data_t& lphi = a_dataPtrsDst[dstAcc.m_dataID][dstAcc.m_offset];
      //always overwriting here.
      lphi = 0;
      for(uint64_t isten = 0; isten < stenSize; isten++)
      {
        const pair_t<data_t>& entry = a_ebstencil[startInd + isten];
        const EBDataLoc& srcAcc = entry.first;
        const data_t   & weight = entry.second;
        const data_t& phi = a_dataPtrsSrc[srcAcc.m_dataID][srcAcc.m_offset];
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
            const data_t                              & a_scale)

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
    aggStencilIndexer<<<m_stride, m_blocks>>>
      (begin, end, dataPtrsSrc, dataPtrsDst, d_ebstencil, d_stensizes, d_stenstart, d_dstaccess);
  }
#endif
}
///////////
#endif
