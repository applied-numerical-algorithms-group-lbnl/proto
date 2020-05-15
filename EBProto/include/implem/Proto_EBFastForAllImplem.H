
#ifndef _Proto_EBFastForAllImplem_H_
#define _Proto_EBFastForAllImplem_H_

namespace Proto
{
///
  template<typename T>
  inline T
  fastGetVar(unsigned int ivec,  T a_s)
  {
    return a_s;
  }

///
  template<CENTERING cent, typename data_t, unsigned int ncomp>
  inline Var<data_t, ncomp>
  fastGetVar(unsigned int a_ivec,
             shared_ptr< vector< EBIrregStruct<cent, data_t, ncomp> >  >a_dst)
  {
    auto dst = *a_dst;
    Var<data_t, ncomp> retval;
    const EBIrregStruct<cent, data_t, ncomp>& ugly = dst[a_ivec];
    for(int icomp = 0; icomp < ncomp; icomp++)
    {
      retval.m_ptrs[icomp] = ugly.m_startPtr + ugly.m_offset + (ugly.m_varsize*icomp);
    }

    return retval;
  }
  
  template<CENTERING cent, typename data_t,unsigned int ncomp,  typename Func, typename... Srcs>
  inline void
  hostFastVectorFunc(const Func& a_F, shared_ptr< vector< EBIrregStruct<cent, data_t, ncomp> > > a_dst, Srcs... a_srcs)
  {
    for(unsigned int ivec = 0; ivec < a_dst.size(); ivec++)
    {
      a_F(fastGetVar(ivec, a_dst), (getVar(ivec, a_srcs))...);
    }
       
  }


///going into this srcs are vector<EBIrregStruct<cent, data_t, ncomp> > and other stuff
  template<CENTERING cent, typename data_t,unsigned int ncomp,  typename Func, typename... Srcs>
  inline void
  hostFastVectorFunc_i(const Func& a_F, shared_ptr<vector< EBIrregStruct<cent, data_t, ncomp> > > a_dst, Srcs... a_srcs)
  {
    auto dst = *a_dst;
    for(unsigned int ivec = 0; ivec < dst.size(); ivec++)
    {
      Point pt = dst[ivec].m_index;
      a_F(pt.m_tuple, fastGetVar(ivec, a_dst), (fastGetVar(ivec, a_srcs))...);
    }
       
  }
  template <typename T>
  inline T
  getFastEBIrregStruct(T& a_T)
  {
    return a_T;
  }

  template <CENTERING cent, typename  data_t, unsigned int ncomp>
  inline std::shared_ptr< vector< EBIrregStruct<cent, data_t, ncomp> > >
  getFastEBIrregStruct(IrregData<cent, data_t, ncomp>& a_s )
  {
    
    return a_s.getUglyStructPtr();
  }

///going into this srcs are IrregDatas and other stuff
  template<CENTERING cent, typename  data_t, unsigned int ncomp, typename Func, typename... Srcs>
  inline void
  hostEBFastForAllIrreg(const Func& a_F, 
                        IrregData<cent, data_t, ncomp>& a_dst,
                        Srcs&...  a_srcs)
  {
    //this uses all the indicies of the input data
    hostFastVectorFunc(a_F, getFastEBIrregStruct(a_dst), (getFastEBIrregStruct(a_srcs))...);
  }


///going into this srcs are IrregDatas and other stuff
  template<CENTERING cent, typename  data_t, unsigned int ncomp, typename Func, typename... Srcs>
  inline void
  hostEBFastForAllIrreg_i(const Func& a_F,
                          IrregData<cent, data_t, ncomp>& a_dst,
                          Srcs&...  a_srcs)
  {
    //this uses all the indicies of the input data
    hostFastVectorFunc_i(a_F, getFastEBIrregStruct(a_dst), (getFastEBIrregStruct(a_srcs))...);
  }


///going into this srcs are EBBoxDatas and other stuff
  template<typename Func, typename... Srcs>
  inline void
  hostEBFastforall(const Func & a_F,  Box a_box, Srcs&... a_srcs)
  {
//call regular forall
    {
      //getbox data livews in Proto_EBForAllImplem
      PR_TIME("regular_forall");
      forallInPlaceBase(a_F, a_box, (getBoxData(a_srcs))...);
    }  
//do the same thing for the irregular data 
    //get irreg data livews in Proto_EBForAllImplem
    {
      PR_TIME("irreg_forall");
      hostEBFastForAllIrreg(a_F, (getIrregData(a_srcs))...);
    }
  }

///going into this srcs are EBBoxDatas and other stuff
  template<typename Func, typename... Srcs>
  inline void
  hostEBFastforall_i(const Func & a_F,  Box a_box, Srcs&... a_srcs)
  {
    {
      //getbox data livews in Proto_EBForAllImplem
      PR_TIME("regular_forall_i");
      forallInPlace_i(a_F, a_box, (getBoxData(a_srcs))...);
    }  
    {
      //get irreg data livews in Proto_EBForAllImplem
      PR_TIME("irreg_forall_i");
      hostEBFastForAllIrreg_i(a_F, (getIrregData(a_srcs))...);
    }
  }
#ifdef PROTO_CUDA
//going into this srcs are EBIrregStruct<cent, data_t, ncomp>* and other stuff
  template<CENTERING cent, typename data_t,unsigned int ncomp,  typename Func, typename... Srcs>
  inline void
  cudaFastVectorFunc(const Func& a_F, unsigned int a_Nvec, 
                 EBIrregStruct<cent, data_t, ncomp> * a_dst,Srcs... a_srcs)
  {
    cudaStream_t curstream = DisjointBoxLayout::getCurrentStream();
    const int N = a_Nvec;
    unsigned int stride = a_Nvec;
    unsigned int blocks = 1;
    size_t smem = 0;
    vec_indexer<<<blocks, stride, smem, curstream>>>
      (0, N, mapper(a_F), a_dst, a_srcs...);

  }


//going into this srcs are EBIrregStruct<cent, data_t, ncomp>* and other stuff
  template<CENTERING cent, typename data_t,unsigned int ncomp,  typename Func, typename... Srcs>
  inline void
  cudaFastVectorFunc_i(const Func& a_F, unsigned int a_Nvec, 
                   EBIrregStruct<cent, data_t, ncomp> * a_dst,Srcs... a_srcs)
  {
    cudaStream_t curstream = DisjointBoxLayout::getCurrentStream();
    const int N = a_Nvec;
    unsigned int stride = a_Nvec;
    unsigned int blocks = 1;
    size_t smem = 0;
    vec_indexer_i<<<blocks, stride, smem, curstream>>>
      (0, N, mapper(a_F), a_dst, a_srcs...);

  }

  template <typename T>
  inline T
  cudagetFastEBIrregStruct(T& a_T)
  {
    return a_T;
  }
//
  template <CENTERING cent, typename  data_t, unsigned int ncomp>
  inline  EBIrregStruct<cent, data_t, ncomp>*
  cudagetFastEBIrregStruct(IrregData<cent, data_t, ncomp>& a_s )
  {
    EBIrregStruct<cent, data_t, ncomp>* retval = a_s.getUglyStructDevicePtr();
    return retval;
  }

//going into this srcs are IrregDatas and other stuff
  template<CENTERING cent, typename  data_t, unsigned int ncomp, typename Func, typename... Srcs>
  inline void
  cudaEBFastForAllIrreg(const Func& a_F, 
                        IrregData<cent, data_t, ncomp>& a_dst,
                        Srcs&...  a_srcs)
  {
    //uses all the incices of a_dst
    int vecsize = a_dst.vecsize();
    cudaVectorFunc(a_F, vecsize, cudagetFastEBIrregStruct(a_dst), 
                   (cudagetFastEBIrregStruct(a_srcs))...);
  }


///going into this srcs are IrregDatas and other stuff
  template<CENTERING cent, typename  data_t, unsigned int ncomp, typename Func, typename... Srcs>
  inline void
  cudaEBFastForAllIrreg_i(const Func& a_F,
                          IrregData<cent, data_t, ncomp>& a_dst,
                          Srcs&...  a_srcs)
  {
    //uses all the incices of a_dst
    int vecsize = a_dst.vecsize();
    cudaFastVectorFunc_i(a_F, vecsize, cudagetFastEBIrregStruct(a_dst), 
                         (cudagetFastEBIrregStruct(a_srcs))...);

  }

  ///
  template<typename Func, typename... Srcs>
  inline void
  cudaEBFastforall(const Func & a_F,  Box a_box, Srcs&... a_srcs)
  {
//call regular forall
    forallInPlaceBase(a_F, a_box, (getBoxData(a_srcs))...);
  
//do the same thing for the irregular data
    cudaFastEBForAllIrreg(a_F, getIrregData(a_srcs)...);
  }


  template<typename Func, typename... Srcs>
  inline void
  cudaEBFastforall_i(const Func & a_F,  Box a_box, Srcs&... a_srcs)
  {
//call regular forall
    forallInPlace_i(a_F, a_box, (getBoxData(a_srcs))...);
  
//do the same thing for the irregular data
    cudaEBFastForAllIrreg_i(a_F, getIrregData(a_srcs)...);
  }
#endif

///version that does not send the point to the function
  template<typename Func, typename... Srcs>
  inline void ebFastforallInPlace(unsigned long long int a_num_flops_point,
                                  const char*            a_timername,
                                  const Func & a_F,  Box a_box, Srcs&... a_srcs)
  {

    PR_TIME(a_timername);

//  printf("in ebforall function pointer = %p\n", &a_F);
    unsigned long long int boxfloops = a_num_flops_point*a_box.size();

#ifdef PROTO_CUDA
    cudaEBFastforall(a_F, a_box, a_srcs...);
    cudaDeviceSynchronize();
#else
    hostEBFastforall(a_F, a_box, a_srcs...);
#endif
    PR_FLOPS(boxfloops);
  }


///version that sends the point to the function
  template<typename Func, typename... Srcs>
  inline void ebFastforallInPlace_i(unsigned long long int a_num_flops_point,
                                    const char*            a_timername,
                                    const Func & a_F,  Box a_box, Srcs&... a_srcs)
  {
    PR_TIME(a_timername);

    unsigned long long int boxfloops = a_num_flops_point*a_box.size();

#ifdef PROTO_CUDA
    cudaEBFastforall_i(a_F, a_box, a_srcs...);
#else
    hostEBFastforall_i(a_F, a_box, a_srcs...);
#endif


    PR_FLOPS(boxfloops);
  }

  ///version that only touches cut cells
  template<typename Func, typename... Srcs>
  inline void ebFastforallIrreg(const char*            a_timername,
                                const Func & a_F,  Box a_box, Srcs&... a_srcs)

  {
    PR_TIME(a_timername);
    //the flops thing here is tricky
#ifdef PROTO_CUDA
    cudaEBFastForAllIrreg(a_F,        getIrregData(a_srcs)...);
#else
    hostEBFastForAllIrreg(a_F,  (getIrregData(a_srcs))...);
#endif
  }



  ///version that only touches cut cells
  template<typename Func, typename... Srcs>
  inline void ebFastforallIrreg_i(const char*            a_timername,
                                  const Func & a_F,     Srcs&... a_srcs)

  {
    PR_TIME(a_timername);
    //the flops thing here is tricky
#ifdef PROTO_CUDA
    cudaEBFastForAllIrreg_i(a_F,  getIrregData(a_srcs)...);
#else
    hostEBFastForAllIrreg_i(a_F,  (getIrregData(a_srcs))...);
#endif
  }
}

#endif

