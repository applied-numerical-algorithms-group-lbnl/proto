
#ifndef _Proto_EBBoxDATAImplem_H_
#define _Proto_EBBoxDATAImplem_H_

namespace Proto
{
  template <CENTERING cent, typename T>
  inline T
  getFastEBIrregStruct(T& a_T)
  {
    return a_T;
  }

  template <CENTERING cent, typename  data_t, unsigned int ncomp>
  inline vector< EBIrregStruct<cent, data_t, ncomp> >*
  getFastEBIrregStruct(IrregData<cent, data_t, ncomp>& a_s )
  {
    
    return a_s.getUglyStructPtr();
  }

///going into this srcs are IrregDatas and other stuff
  template<CENTERING cent, typename  data_t, unsigned int ncomp, typename Func, typename... Srcs>
  inline void
  hostEBFastForAllIrreg(const Func& a_F, const Box& a_box,
                        IrregData<cent, data_t, ncomp>& a_dst,
                        Srcs&...  a_srcs)
  {
    //this uses all the indicies of the input data
    hostVectorFunc(a_F, getFastEBIrregStruct(dstvofs, a_dst), (getFastEBIrregStruct(dstvofs, a_srcs))...);
  }


///going into this srcs are IrregDatas and other stuff
  template<CENTERING cent, typename  data_t, unsigned int ncomp, typename Func, typename... Srcs>
  inline void
  hostEBFastForAllIrreg_i(const Func& a_F, const Box& a_box,
                          IrregData<cent, data_t, ncomp>& a_dst,
                          Srcs&...  a_srcs)
  {
    //this uses all the indicies of the input data
    hostVectorFunc_i(a_F, getFastEBIrregStruct(a_dst), (getFastEBIrregStruct(a_srcs))...);
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
      hostEBFastForAllIrreg_i(a_F, getIrregData(a_srcs)...);
    }
  }


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

