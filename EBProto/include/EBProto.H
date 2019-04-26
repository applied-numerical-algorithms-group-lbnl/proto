#pragma once


#include "Proto.H"
#include "Proto_GeometryService.H"
#include "Proto_PointSet.H"
#include "Proto_HostBoxData.H"
#include "Proto_SimpleImplicitFunctions.H"

namespace Proto
{
  enum CENTERING {CELL,BOUNDARY,XFACE,YFACE,ZFACE};

/// Return the index into the arc vector
  inline int ebp_index(int a_idir, Side::LoHiSide a_side) 
  {
    return a_idir + a_sd*DIM;
  }
}
