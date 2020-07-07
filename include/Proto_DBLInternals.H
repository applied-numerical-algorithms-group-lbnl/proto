#ifndef _PROTO_DBLInternals_H_
#define _PROTO_DBLInternals_H_

#include "Proto_Point.H"
#include "Proto_Box.H"
#include "implem/Proto_Stack.H"
#include "Proto_ProblemDomain.H"
#include "Proto_Morton.H"
#include "Proto_SPMD.H"
#include <cstdlib> //for size_t
#include <iostream>
#include <stack>
#include <memory>
#include <array>
#include <unordered_map>
namespace Proto 
{
  typedef std::unordered_map<Point, int, Proto::MortonOrder > BoxMap;
  struct DBLInternals
    {
      BoxMap                    boxMap;
      // BoxMap stores an int pointing to a location in allBoxes, indexed by Points.
      vector<pair<Point, int> > allBoxes;
      // allBoxes stores a Point corresponding to the box, and the processor number. Locations
      // for any given processor are stored in a contiguous set.
      vector<int>               startProc;
      // Stores the starting location in allBoxes for each processor.
      ProblemDomain             coarsenedDomain;
      // m_problemDomain coarsened by m_boxSize. Remains invariant under coarsening or refinement
      // so we can put it here.
      
    };
}
#endif
