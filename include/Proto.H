#pragma once

#ifdef PROTO_CUDA

#define CUDA_DECORATION __host__ __device__

#else

#define CUDA_DECORATION /*__host__ __device__ */ 

#endif


#include <cstdlib>
#include "Proto_RK4.H"
#include "Proto_PAssert.H"
#include "Proto_PowerItoI.H"
#include "Proto_Point.H"
#include "Proto_Box.H"
#include "Proto_BoxData.H"
#include "Proto_Stencil.H"
#include "Proto_LevelData.H"
//#include "WriteBoxData.H"
