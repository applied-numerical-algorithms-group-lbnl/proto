#ifndef __PROTO_MMB_HEADER__
#define __PROTO_MMB_HEADER__

#ifndef PR_MMB
#define PR_MMB
#endif

//#include "Proto.H"

#include "mmb/Proto_MBGraph.H"
#include "mmb/Proto_MBProblemDomain.H"
#include "mmb/Proto_MBBoxPartition.H"
#include "mmb/Proto_MBDisjointBoxLayout.H"
#include "mmb/Proto_MBLevelBoxData.H"
#include "mmb/Proto_MBLevelMap.H"
#include "mmb/Proto_MBLevelOp.H"
#ifdef PR_OPS 
#include "mmb/Proto_MBInterpOp.H"
#endif
#ifdef PR_AMR
#include "ProtoAMR.H"
#include "mmb/Proto_MBAMRBoxPartition.H"
#include "mmb/Proto_MBAMRGrid.H"
#endif
#endif
