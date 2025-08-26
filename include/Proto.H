#ifndef __PROTO_HEADER__
#define __PROTO_HEADER__

// UTILITY
#include "base/Proto_Timer.H"
#include "base/Proto_PAssert.H"
#include "base/Proto_SPMD.H"
#include "base/Proto_Parstream.H" //SPMD
#include "base/Proto_MayDay.H" //parstream
#include "base/Proto_MemType.H"
#include "base/Proto_Memory.H" //MemType

#include "base/Proto_CPU.H" //MemInfo
#include "base/Proto_GPU.H" //MemInfo
#include "base/Proto_MemInfo.H" //MemType
#include "base/Proto_Macros.H" //CPU, GPU, MemType
#include "base/Proto_Stack.H" // MemType, Macros, Memory

#include "base/implem/Proto_TimerImplem.H"
#include "base/implem/Proto_ParstreamImplem.H"
#include "base/implem/Proto_MayDayImplem.H"
#include "base/implem/Proto_MemoryImplem.H"

// SINGLE PATCH OBJECTS
#ifdef PR_OPS
#include "ProtoOps.H"
#endif
#include "base/Proto_Accel.H"
#include "base/Proto_Math.H"
#include "base/Proto_Face.H"
#include "base/Proto_Centering.H"
#include "base/Proto_CInterval.H"
#include "base/Proto_Point.H" // [Accel], Face
#include "base/Proto_Array.H" // [Accel], Point
#include "base/Proto_Box.H" // Face, Point, Centering
#include "base/Proto_CoordPermutation.H" //Point, Array, Box
#include "base/Proto_Reduction.H" //SPMD, Memory, Macros, [Accel]
#include "base/Proto_BoxData.H" //MemType, Memory, Stack, Array, Box, CoordPermutation, CInterval, Reduction
#include "base/Proto_Stencil.H" //MemType, Point, BoxData, Math, [Matrix]
#include "base/Proto_Forall.H" //MemType, Point, Box, BoxData

#include "base/implem/Proto_PointImplem.H"
#include "base/implem/Proto_ArrayImplem.H"
#include "base/implem/Proto_BoxImplem.H"
#include "base/implem/Proto_CoordPermutationImplem.H"
#include "base/implem/Proto_CIntervalImplem.H"
#include "base/implem/Proto_ReductionImplem.H"
#include "base/implem/Proto_BoxDataImplem.H"
#include "base/implem/Proto_StencilImplem.H"
#include "base/implem/Proto_StencilDefs.H"
#include "base/implem/Proto_ForallImplem.H"




#include "base/Proto_DistributedSet.H"


#ifdef PR_HDF5
#include "base/Proto_HDF5.H"
#endif


#include "base/Proto_InterpStencil.H"
#include "base/Proto_BoundaryCondition.H"
#include "base/Proto_BoxOp.H"
#include "base/Proto_ProblemDomain.H"
#include "base/Proto_BoxPartition.H"
#include "base/Proto_DataIndex.H"
#include "base/Proto_DataIterator.H"
#include "base/Proto_DisjointBoxLayout.H"
#include "base/Proto_LevelBoxData.H"

#include "base/Proto_CompositeOp.H"
#include "base/Proto_Operator.H"
#include "base/Proto_LevelBC.H"
#include "base/Proto_LevelOp.H"


#ifdef PR_AMR
#include "ProtoAMR.H"
#endif
#ifdef PR_MMB
#include "ProtoMMB.H"
#endif

#include "transition/Proto_Transition_Utilities.H"

#include "utils/Proto_LevelMultigridSolver.H"
#include "utils/Proto_LevelRK4.H"

#endif
