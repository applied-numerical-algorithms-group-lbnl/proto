#ifndef __PROTO_CUDA_HEADER__
#define __PROTO_CUDA_HEADER__


#ifdef PROTO_CUDA

#define CUDA_DECORATION __host__ __device__

#else

#define CUDA_DECORATION /*__host__ __device__ */ 

#endif


//#include "WriteBoxData.H"
#endif
