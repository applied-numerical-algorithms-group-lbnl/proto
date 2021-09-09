#ifndef __PROTO_HEADER__
#define __PROTO_HEADER__

#include <cstdlib>
#include "Proto_parstream.H"
#include "Proto_macros.H"
#include "Proto_SPMD.H"
#include "Proto_MemInfo.H"
#ifdef PROTO_CUDA
#include "Proto_cuda.H"
#endif
#include "Proto_Point.H"
#include "Proto_Box.H"
#include "Proto_RK4.H"
#include "Proto_PAssert.H"
#include "Proto_PowerItoI.H"
#include "Proto_BoxData.H"
#include "Proto_Stencil.H"
#include "Proto_DisjointBoxLayout.H"
#include "Proto_LevelBoxData.H"
#include "Proto_Reduction.H"
#include "Proto_AMRGrid.H"
#include "Proto_AMRData.H"
#include "Proto_FinitePointSet.H"
#ifdef PR_HDF5
#include "Proto_HDF5.H"
#endif
#include "Proto_MayDay.H"
#include "implem/Proto_FusedStencil.H"
#endif
