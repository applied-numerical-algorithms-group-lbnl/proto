#ifndef _PROTO_DisjointBoxLayout_H_
#define _PROTO_DisjointBoxLayout_H_

#include "Proto_Point.H"
#include "Proto_Box.H"
#include "implem/Proto_Stack.H"
#include "Proto_ProblemDomain.H"
#include "Proto_DBLInternals.H"
#include "Proto_DataIndex.H"
#include "Proto_Morton.H"
#include "Proto_SPMD.H"
#include <cstdlib> //for size_t
#include <iostream>
#include <stack>
#include <memory>
#include <array>
#include <unordered_map>

using namespace std;
using std::shared_ptr;
using std::array;
using std::vector;
using std::stack;
 
namespace Proto 
{
  /**
   Simplified version Chombo's DBL for fixed-box sizes.
  */
  // forward declaration of iterators.
  class DisjointBoxLayout
  {
    
    friend class DataIndex;
  private:
    ProblemDomain             m_problemDomain;
    Point                     m_boxSize;    
    shared_ptr<DBLInternals>  m_internals;
    // This makes the proliferation of copies not be a performance barrier, as well as 
    // making the data pointed at by m_internals reuseable under coarsening or refinement.
    // If making m_problemDomain and m_boxSize objects becomes a performance issue, can wrap 
    // them in a shared pointer as well.
  public:
    ///
    inline DisjointBoxLayout();

    /// Constructor for an entire rectangle being distributed over processors.
    /** Uses Morton ordering to distribute boxes to processors. */
     
    inline DisjointBoxLayout(const ProblemDomain   & a_problemDomain, 
                             const Point           & a_boxSize);

    /// Constructor for a collection of Boxes distributed over processors.
    /** Uses Morton-ordering to distribute boxes to processors. */
    
    inline DisjointBoxLayout(const ProblemDomain   & a_problemDomain,
                             const vector<Point>   & a_coarsenedPatches,
                             const Point           & a_boxSize);
 
    ///
    inline DisjointBoxLayout(const DisjointBoxLayout& a_input)
    {
      m_internals = a_input.m_internals;
      m_problemDomain = a_input.m_problemDomain;
      m_boxSize = a_input.m_boxSize;
    };
    ///
    inline void 
    define(const ProblemDomain   & a_problemDomain, 
           const Point           & a_boxSize);
    ///
    inline void
    define(const ProblemDomain   & a_problemDomain,
           const vector<Point>   & a_coarsenedPatches,
           const Point           & a_boxSize);
    ///
    inline DisjointBoxLayout& operator=(const DisjointBoxLayout& a_input)
    {
      m_internals = a_input.m_internals;
      m_problemDomain = a_input.m_problemDomain;
      m_boxSize = a_input.m_boxSize;
      return *this;
    };
    ///
    inline bool operator==(const DisjointBoxLayout& a_input) const
    {
      return ((a_input.m_internals == m_internals)&&
              (a_input.m_problemDomain == m_problemDomain)&&
              (a_input.m_boxSize == m_boxSize));
    };
    /*
    inline void dump() const
    {
      using std::cout;
      using std::endl;
      cout << "dbl domain   = " << m_problemDomain  << endl;
      cout << "dbl coarse domain  = " << m_internals->coarsenedDomain << endl;
      cout << "dbl box size = " << m_boxSize << endl;
      // need to use iterator over m_procmap.
      cout << "dbl boxes: "   << endl;
      for(int ibox = 0; ibox  < size(); ibox++)
      {
        cout << ibox << " " << (*this)[ibox]      << endl;
        } */
    /// Indexing with a DataIndex - returns a Point corresponding to the coarsened patch.
    inline Point operator[](const DataIndex& a_index) const;

    /// Get procID
    inline int procID(const DataIndex& a_index) const;

    /// Finding the DataIndex corresponding to the Point in the bitmap defining the DBL.
    /** This includes Points outside the ProblemDomain Box, provided they are in periodic 
        images of that Box. */
    inline DataIndex index(Point& a_coarsenedPt) const;
    
    ///number of boxes in dbl. 
    inline unsigned int size() const;

    ///number of boxes in dbl on myProc. 
    inline unsigned int mySize() const;

    /// Indexing into the Boxes in the DBL.
    inline DataIndex index(unsigned int a_indexInt) const;
    
    /// Indexing into the Boxes on myProc.
    inline DataIndex myIndex(unsigned int a_myIndexInt) const;

    /// Indexing into the Boxes on myProc.
    inline int myIntIndex(const DataIndex& a_di) const;

    /// Get box size.
    inline Point boxSize() const;

    /// Get Problem Domain.
    inline ProblemDomain problemDomain() const;

    /// Find function: if the Point is in the DBL, return the DataIndex for that
    /// Box; otherwise, return an invalid DataIndex.
    inline DataIndex find(const Point& a_pt) const;
    
    /// End function for constructing iterators.
    inline const DataIndex& end() const;
    
    /// 
    inline bool coarsenable(const Point& a_refrat) const;

    /// Returns DBL coarsened by a_refrat.
    inline DisjointBoxLayout coarsen(const Point& a_refrat) const;

    /// Returns DBL refined by a_refrat. 
    inline DisjointBoxLayout refine(const Point& a_refrat) const;
  };
#include "implem/Proto_DisjointBoxLayoutImplem.H"
} //end namespace Proto

#endif
