#ifndef __PROTO_HEADER__
#define __PROTO_HEADER__



#include <cstdlib>
#ifdef PROTO_CUDA
#include "Proto_gpu.H"
#include "Proto_cuda.H"
#else
#include "Proto_cpu.H"
#endif
#include "Proto_Point.H"
#include "Proto_Box.H"
#include "Proto_RK4.H"
#include "Proto_PAssert.H"
#include "Proto_PowerItoI.H"
#include "Proto_BoxData.H"
#include "Proto_Stencil.H"
#include "implem/Proto_DisjointBoxLayout.H"
#include "implem/Proto_LevelData.H"
#include "Proto_Reduction.H"
//#include "WriteBoxData.H"
#include "implem/Proto_FusedStencil.H"
#endif
