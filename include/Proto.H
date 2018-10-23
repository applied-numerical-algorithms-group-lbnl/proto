#pragma once

#ifdef PROTO_CUDA

#define CUDA_DECORATION __host__ __device__

#else

#define CUDA_DECORATION /*__host__ __device__ */ 

#endif


#include <cstdlib>
#include "RK4.H"
#include "PAssert.H"
#include "PowerItoI.H"
#include "Point.H"
#include "Bx.H"
#include "BoxData.H"
#include "Stencil.H"
//#include "WriteBoxData.H"
