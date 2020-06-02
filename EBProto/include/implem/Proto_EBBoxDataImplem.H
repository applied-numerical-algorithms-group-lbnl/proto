
#ifndef _Proto_EBBoxDATAImplem_H_
#define _Proto_EBBoxDATAImplem_H_

namespace Proto
{

  template<typename T>
  inline T&
  getIrregData(T& a_s)
  {
    return a_s;
  }

  template <CENTERING cent, typename  data_t, unsigned int ncomp>
  inline IrregData<cent, data_t, ncomp>&
  getIrregData(EBBoxData<cent, data_t, ncomp>& a_s)
  {
    return a_s.getIrregData();
  }



  template <CENTERING cent, typename  data_t, unsigned int ncomp>
  inline IrregData<cent, data_t, ncomp>&
  getIrregData(const EBBoxData<cent, data_t, ncomp>& a_s)
  {
    return (IrregData<cent, data_t, ncomp>&)(a_s.getIrregData());
  }


  template<CENTERING cent,  typename data_t, unsigned int ncomp>
  struct
  EBIrregStruct
  {
    data_t*        m_startPtr;
    unsigned int   m_varsize;
    unsigned int   m_offset;
    Point          m_index;
  };

  template <CENTERING cent, typename T>
  inline T
  getEBIrregStruct(const vector<EBIndex<cent> >& a_indices,
                T& a_T)
  {
    return a_T;
  }

  template <CENTERING cent, typename  data_t, unsigned int ncomp>
  inline vector< EBIrregStruct<cent, data_t, ncomp> >
  getEBIrregStruct(const vector<EBIndex<cent> >& a_indices,
                   IrregData<cent, data_t, ncomp>& a_s )
  {
//    data_t*      debPtr  = a_s.data();
//  printf("irreg data ptr = %p\n", debPtr);
    vector< EBIrregStruct<cent, data_t, ncomp> > retval;
    for(int ivec = 0; ivec < a_indices.size(); ivec++)
    {
////begin debug
//      Point pt = a_indices[ivec].m_pt;
//      if((pt[0]==17) && pt[1]==18)
//      {
//        std::cout << "pt = " << pt << ", ivec = " << ivec << std::endl;
//      }
////end debug
      EBIrregStruct<cent, data_t, ncomp>  vecval;
      vecval.m_startPtr = a_s.data();
      vecval.m_varsize  = a_s.vecsize();
      vecval.m_offset   = a_s.index(a_indices[ivec], 0);
      vecval.m_index    = a_indices[ivec].m_pt;
      retval.push_back(vecval);
    }

//  printf("host return data ptr = %p\n", retval.data());
    return retval;
  }


  template<typename T>
  inline T&
  getBoxData(T& a_s)
  {
    return a_s;
  }

  template <CENTERING cent, typename  data_t, unsigned int ncomp>
  inline BoxData<data_t, ncomp>&
  getBoxData(EBBoxData<cent, data_t, ncomp>& a_s)
  {
    return a_s.getRegData();
  }


  template <CENTERING cent, typename  data_t, unsigned int ncomp>
  inline const BoxData<data_t, ncomp>&
  getBoxData(const EBBoxData<cent, data_t, ncomp>& a_s)
  {
    return a_s.getRegData();
  }


#ifdef PROTO_CUDA

  template<typename T>
  __device__ __host__
  inline T
  cudaGetVar(unsigned int ivec,  T a_s)
  {
    return a_s;
  }


  template<CENTERING cent, typename data_t, unsigned int ncomp>
  __device__ 
  inline Var<data_t, ncomp>
  cudaGetVar(unsigned int a_ivec,
             EBIrregStruct<cent, data_t, ncomp>* a_dst)
  {
    Var<data_t, ncomp> retval;

    const EBIrregStruct<cent, data_t, ncomp>*  rawptr = a_dst;
    const EBIrregStruct<cent, data_t, ncomp>&  ugly   = rawptr[a_ivec];
//  printf("cudaGetVar rawptr   = %p \n", rawptr);
//  printf("cudaGetVar startptr = %p \n", ugly.m_startPtr);
    for(int icomp = 0; icomp < ncomp; icomp++)
    {
      retval.m_ptrs[icomp] = ugly.m_startPtr + ugly.m_offset + (ugly.m_varsize*icomp);
    }
    return retval;
  }




//going into this srcs are thrust_device_pointer<uglystruct> and other stuff
  template<CENTERING cent, typename data_t,unsigned int ncomp,  typename Func, typename... Srcs>
  __global__
  void
  vec_indexer(unsigned int a_begin, unsigned int a_end,Func a_body, 
              EBIrregStruct<cent, data_t, ncomp>*  a_dst, Srcs... a_srcs)
  {
    int idx = threadIdx.x + blockIdx.x*blockDim.x;
    if (idx >= a_begin && idx < a_end)
    {
      a_body(cudaGetVar(idx, a_dst), cudaGetVar(idx, a_srcs)...);
    }
  }


  template<CENTERING cent, typename data_t,unsigned int ncomp,  typename Func, typename... Srcs>
  __global__
   void
  vec_indexer_i(unsigned int a_begin, unsigned int a_end,Func a_body, 
                EBIrregStruct<cent, data_t, ncomp>*  a_dst, Srcs... a_srcs)
  {
    int idx = threadIdx.x + blockIdx.x*blockDim.x;
    if (idx >= a_begin && idx < a_end)
    {
      a_body(a_dst[idx].m_index.m_tuple, cudaGetVar(idx, a_dst), cudaGetVar(idx, a_srcs)...);
    }
  }
//going into this srcs are uglystruct* and other stuff
  template <typename T>
  inline int
  cleanUpPtrMem(T a_T)
  {
    return 0;
  }
//going into this srcs are uglystruct* and other stuff
  template <CENTERING cent, typename data_t,unsigned int ncomp>
  inline int
  cleanUpPtrMem(EBIrregStruct<cent, data_t, ncomp> * a_ptr)
  {
    protoFree(a_ptr);
    return 0;
  }

//
  template<typename... Srcs>
  inline void
  emptyFunc(Srcs... a_srcs)
  {
  }

//going into this srcs are uglystruct* and other stuff
  template<CENTERING cent, typename data_t,unsigned int ncomp,  typename Func, typename... Srcs>
  inline void
  cudaVectorFunc(const Func& a_F, unsigned int a_Nvec, 
                 EBIrregStruct<cent, data_t, ncomp> * a_dst,Srcs... a_srcs)
  {
//  printf("cudavecf: dst  = %p\n", a_dst);
    //printf("cudavecf: src  = %p\n", a_firstsrc);
    protoStream_t curstream = DisjointBoxLayout::getCurrentStream();
    const int N = a_Nvec;
    unsigned int stride = a_Nvec;
    unsigned int blocks = 1;
    size_t smem = 0;
    protoLaunchKernelMemAsync(vec_indexer, blocks, stride, smem, curstream,
      0, N, mapper(a_F), a_dst, a_srcs...);

    //there is a cudaMalloc that happens above so we have to delete
    emptyFunc(cleanUpPtrMem(a_dst ), (cleanUpPtrMem(a_srcs))...); 
  }


//going into this srcs are uglystruct* and other stuff
  template<CENTERING cent, typename data_t,unsigned int ncomp,  typename Func, typename... Srcs>
  inline void
  cudaVectorFunc_i(const Func& a_F, unsigned int a_Nvec, 
                   EBIrregStruct<cent, data_t, ncomp> * a_dst,Srcs... a_srcs)
  {
//  printf("cudavecf: dst  = %p\n", a_dst);
    //printf("cudavecf: src  = %p\n", a_firstsrc);
    protoStream_t curstream = DisjointBoxLayout::getCurrentStream();
    const int N = a_Nvec;
    unsigned int stride = a_Nvec;
    unsigned int blocks = 1;
    size_t smem = 0;
    protoLaunchKernelMemAsync(vec_indexer_i, blocks, stride, smem, curstream,
      0, N, mapper(a_F), a_dst, a_srcs...);

    //there is a cudaMalloc that happens above so we have to delete
    emptyFunc(cleanUpPtrMem(a_dst ), (cleanUpPtrMem(a_srcs))...); 
  }


  template <CENTERING cent, typename T>
  inline T
  cudagetEBIrregStruct(const vector<EBIndex<cent> >& a_indices,
                    T& a_T)
  {
    return a_T;
  }
//
  template <CENTERING cent, typename  data_t, unsigned int ncomp>
  inline  EBIrregStruct<cent, data_t, ncomp>*
  cudagetEBIrregStruct(const vector<EBIndex<cent> >& a_indices,
                    IrregData<cent, data_t, ncomp>& a_s )
  {
    vector< EBIrregStruct<cent, data_t, ncomp> > hostvec = getEBIrregStruct(a_indices, a_s);

    size_t memsize = hostvec.size()*sizeof(EBIrregStruct<cent, data_t, ncomp>);
    EBIrregStruct<cent, data_t, ncomp>* retval;
    protoMalloc(&retval, memsize);
    protoMemcpy(retval, hostvec.data(), memsize, protoMemcpyHostToDevice);

    //this copies from the host to the device
//  printf("cgus: device host vector ptr = %p\n", hostvec.data());
//  printf("cgus: device return data ptr = %p\n", retval);
    return retval;
  }
//going into this srcs are IrregDatas and other stuff
  template<CENTERING cent, typename  data_t, unsigned int ncomp, typename Func, typename... Srcs>
  inline void
  cudaEBForAllIrreg(const Func& a_F, const Box& a_box,
                    IrregData<cent, data_t, ncomp>& a_dst,
                    Srcs&...  a_srcs)
  {
    //indicies into irreg vector that correspond to input box
    const vector<EBIndex<cent> >& dstvofs = a_dst.getIndices(a_box);
    unsigned int vecsize = a_dst.vecsize();
    if(vecsize > 0)
    {
//    printf("cudaebforall: dst  = %p\n", a_dst.data());
      cudaVectorFunc(a_F, vecsize, cudagetEBIrregStruct(dstvofs, a_dst), 
                     (cudagetEBIrregStruct(dstvofs, a_srcs))...);

    }
  }


///going into this srcs are IrregDatas and other stuff
  template<CENTERING cent, typename  data_t, unsigned int ncomp, typename Func, typename... Srcs>
  inline void
  cudaEBForAllIrreg_i(const Func& a_F, const Box& a_box,
                      IrregData<cent, data_t, ncomp>& a_dst,
                      Srcs&...  a_srcs)
  {
    //indicies into irreg vector that correspond to input box
    const vector<EBIndex<cent> >& dstvofs = a_dst.getIndices(a_box);
    unsigned int vecsize = a_dst.vecsize();
    if(vecsize > 0)
    {
//    printf("cudaebforall: dst  = %p\n", a_dst.data());
      cudaVectorFunc_i(a_F, vecsize, cudagetEBIrregStruct(dstvofs, a_dst), 
                       (cudagetEBIrregStruct(dstvofs, a_srcs))...);

    }
  }

  ///
  template<typename Func, typename... Srcs>
  inline void
  cudaEBforall(const Func & a_F,  Box a_box, Srcs&... a_srcs)
  {
//call regular forall
    forallInPlaceBase(a_F, a_box, (getBoxData(a_srcs))...);
  
//do the same thing for the irregular data
    cudaEBForAllIrreg(a_F, a_box, getIrregData(a_srcs)...);
  }


  template<typename Func, typename... Srcs>
  inline void
  cudaEBforall_i(const Func & a_F,  Box a_box, Srcs&... a_srcs)
  {
//call regular forall
    forallInPlace_i(a_F, a_box, (getBoxData(a_srcs))...);
  
//do the same thing for the irregular data
    cudaEBForAllIrreg_i(a_F, a_box, getIrregData(a_srcs)...);
  }

#else
///cpu-only specific functions

///
  template<typename T>
  inline T
  getVar(unsigned int ivec,  T a_s)
  {
    return a_s;
  }

///
  template<CENTERING cent, typename data_t, unsigned int ncomp>
  inline Var<data_t, ncomp>
  getVar(unsigned int a_ivec,
         vector< EBIrregStruct<cent, data_t, ncomp> > a_dst)
  {
    Var<data_t, ncomp> retval;
    const EBIrregStruct<cent, data_t, ncomp> ugly = a_dst[a_ivec];
    for(int icomp = 0; icomp < ncomp; icomp++)
    {
      retval.m_ptrs[icomp] = ugly.m_startPtr + ugly.m_offset + (ugly.m_varsize*icomp);
    }
//begin debug
//    if(a_ivec == 19)
//    {
//      std::cout << "getvar, ivec == 19, values = "  ;
//      for(int icomp = 0; icomp < ncomp; icomp++)
//      {
//        std::cout << *retval.m_ptrs[icomp] << " ";
//      }
//      std::cout << std::endl;
//    }
//end debug

    return retval;
  }

///going into this srcs are vector<uglystruct> and other stuff
  template<CENTERING cent, typename data_t,unsigned int ncomp,  typename Func, typename... Srcs>
  inline void
  hostVectorFunc(const Func& a_F, vector< EBIrregStruct<cent, data_t, ncomp> > a_dst, Srcs... a_srcs)
  {
    for(unsigned int ivec = 0; ivec < a_dst.size(); ivec++)
    {
      a_F(getVar(ivec, a_dst), (getVar(ivec, a_srcs))...);
    }
       
  }


///going into this srcs are vector<uglystruct> and other stuff
  template<CENTERING cent, typename data_t,unsigned int ncomp,  typename Func, typename... Srcs>
  inline void
  hostVectorFunc_i(const Func& a_F, vector< EBIrregStruct<cent, data_t, ncomp> > a_dst, Srcs... a_srcs)
  {
    for(unsigned int ivec = 0; ivec < a_dst.size(); ivec++)
    {
      Point pt = a_dst[ivec].m_index;
      a_F(pt.m_tuple, getVar(ivec, a_dst), (getVar(ivec, a_srcs))...);
    }
       
  }

///going into this srcs are IrregDatas and other stuff
  template<CENTERING cent, typename  data_t, unsigned int ncomp, typename Func, typename... Srcs>
  inline void
  hostEBForAllIrreg(const Func& a_F, const Box& a_box,
                    IrregData<cent, data_t, ncomp>& a_dst,
                    Srcs&...  a_srcs)
  {
//indicies into irreg vector that correspond to input box
    vector<EBIndex<cent> > dstvofs = a_dst.getIndices(a_box);
    hostVectorFunc(a_F, getEBIrregStruct(dstvofs, a_dst), (getEBIrregStruct(dstvofs, a_srcs))...);
  }


///going into this srcs are IrregDatas and other stuff
  template<CENTERING cent, typename  data_t, unsigned int ncomp, typename Func, typename... Srcs>
  inline void
  hostEBForAllIrreg_i(const Func& a_F, const Box& a_box,
                      IrregData<cent, data_t, ncomp>& a_dst,
                      Srcs&...  a_srcs)
  {
//indicies into irreg vector that correspond to input box
    vector<EBIndex<cent> > dstvofs = a_dst.getIndices(a_box);
    hostVectorFunc_i(a_F, getEBIrregStruct(dstvofs, a_dst), (getEBIrregStruct(dstvofs, a_srcs))...);
  }


///going into this srcs are EBBoxDatas and other stuff
  template<typename Func, typename... Srcs>
  inline void
  hostEBforall(const Func & a_F,  Box a_box, Srcs&... a_srcs)
  {
//call regular forall
    {
      PR_TIME("regular_forall");
      forallInPlaceBase(a_F, a_box, (getBoxData(a_srcs))...);
    }  
//do the same thing for the irregular data
    {
      PR_TIME("irreg_forall");
      hostEBForAllIrreg(a_F, a_box, (getIrregData(a_srcs))...);
    }
  }

///going into this srcs are EBBoxDatas and other stuff
  template<typename Func, typename... Srcs>
  inline void
  hostEBforall_i(const Func & a_F,  Box a_box, Srcs&... a_srcs)
  {
    {
      PR_TIME("regular_forall_i");
      forallInPlace_i(a_F, a_box, (getBoxData(a_srcs))...);
    }  
    {
      PR_TIME("irreg_forall_i");
      hostEBForAllIrreg_i(a_F, a_box, getIrregData(a_srcs)...);
    }
  }

#endif

///version that does not send the point to the function
  template<typename Func, typename... Srcs>
  inline void ebforallInPlace(unsigned long long int a_num_flops_point,
                              const char*            a_timername,
                              const Func & a_F,  Box a_box, Srcs&... a_srcs)
  {
    PR_TIME(a_timername);

//  printf("in ebforall function pointer = %p\n", &a_F);
    unsigned long long int boxfloops = a_num_flops_point*a_box.size();

#ifdef PROTO_CUDA
    cudaEBforall(a_F, a_box, a_srcs...);
    protoDeviceSynchronize();
#else
    hostEBforall(a_F, a_box, a_srcs...);
#endif
    PR_FLOPS(boxfloops);
  }


///version that sends the point to the function
  template<typename Func, typename... Srcs>
  inline void ebforallInPlace_i(unsigned long long int a_num_flops_point,
                                const char*            a_timername,
                                const Func & a_F,  Box a_box, Srcs&... a_srcs)
  {
    PR_TIME(a_timername);

    unsigned long long int boxfloops = a_num_flops_point*a_box.size();

#ifdef PROTO_CUDA
    cudaEBforall_i(a_F, a_box, a_srcs...);
#else
    hostEBforall_i(a_F, a_box, a_srcs...);
#endif


    PR_FLOPS(boxfloops);
  }

  ///version that only touches cut cells
  template<typename Func, typename... Srcs>
  inline void ebforallIrreg(const char*            a_timername,
                            const Func & a_F,  Box a_box, Srcs&... a_srcs)

  {
    PR_TIME(a_timername);
    //the flops thing here is tricky
#ifdef PROTO_CUDA
    cudaEBForAllIrreg(a_F, a_box, getIrregData(a_srcs)...);
#else
    hostEBForAllIrreg(a_F, a_box, (getIrregData(a_srcs))...);
#endif
  }



  ///version that only touches cut cells
  template<typename Func, typename... Srcs>
  inline void ebforallIrreg_i(const char*            a_timername,
                              const Func & a_F,  Box a_box, Srcs&... a_srcs)

  {
    PR_TIME(a_timername);
    //the flops thing here is tricky
#ifdef PROTO_CUDA
    cudaEBForAllIrreg_i(a_F, a_box, getIrregData(a_srcs)...);
#else
    hostEBForAllIrreg_i(a_F, a_box, (getIrregData(a_srcs))...);
#endif
  }
}

#endif

