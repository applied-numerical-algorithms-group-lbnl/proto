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
using namespace Proto;
namespace Proto 
{ 
  /// Simplified version Chombo's DisjointBoxLayout (DBL) for fixed-box sizes.
  class DisjointBoxLayout
  { 
    friend class DataIndex;
  private:
    ProblemDomain             m_problemDomain;
    Point                     m_boxSize;  
    DataIndex                 m_end;     
    shared_ptr<DBLInternals>  m_internals;
    /* 
       The use of  shared_ptr<DBLInternals> makes the proliferation of copies not be a 
       performance barrier, as well as making the data pointed at by m_internals reusable 
       under coarsening or refinement.If making m_problemDomain and m_boxSize objects 
       becomes a performance issue, can wrap them in a shared pointer as well.
    */
  public:
    ///
    inline DisjointBoxLayout();

    /// Constructor for an entire rectangle being distributed over processors.
    /**  a_problemDomain.coarsenable(a_boxSize) must be true.*/
     
    inline DisjointBoxLayout(const ProblemDomain   & a_problemDomain, 
                             const Point           & a_boxSize);

    /// Constructor for a collection of Boxes distributed over processors.
    /** Uses Morton ordering to distribute boxes to processors. 
        \param a_problemDomain   ProblemDomain containing the DBL.
        \param a_coarsenedPatches   contains the Points correcpoinding the patches in the DBL for all the processors.
        \param a_boxSize is the size of a patch. Note that this allows anisotropic patches.

        if ptk = a_coarsenedPatches[k], then Box(ptk,ptk).refine(a_boxSize) is the k-th patch in the DBL.
     */
    
    inline DisjointBoxLayout(const ProblemDomain   & a_problemDomain,
                             const vector<Point>   & a_coarsenedPatches,
                             const Point           & a_boxSize);
 
    /// Copy constructor. Copies shared pointer of m_internals.
    inline DisjointBoxLayout(const DisjointBoxLayout& a_input)
    {
      m_internals = a_input.m_internals;
      m_problemDomain = a_input.m_problemDomain;
      m_boxSize = a_input.m_boxSize;
    };
    /// Define, to permit weak construction. Not sure it is really needed.
    inline void 
    define(const ProblemDomain   & a_problemDomain, 
           const Point           & a_boxSize);
    ///  Define, to permit weak construction. Not sure it is really needed.
    inline void
    define(const ProblemDomain   & a_problemDomain,
           const vector<Point>   & a_coarsenedPatches,
           const Point           & a_boxSize);
    
    /// Assignment operator. Copies shared pointer of m_internals.
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
    /** @name Access operators. 
     */
   
    /// Access to the DBLInternals pointer. Needed to construct DataIndex.
    inline const shared_ptr<DBLInternals> dblInternals() const{return m_internals;};
   
    /// Indexing with a DataIndex - returns a Point corresponding to the coarsened patch.
    inline Point operator[](const DataIndex& a_index) const;
    
    /// Get procid
    inline int procid(const DataIndex& a_index) const;
    
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
    inline Point boxSize() const{return m_boxSize;};
    
    /// Get Problem Domain.
    inline ProblemDomain problemDomain() const{return m_problemDomain;};
    
    /// Find function: if the Point is in the DBL, return the DataIndex for that Box; otherwise, return an invalid DataIndex.
    inline DataIndex find(const Point& a_pt) const;
    ///@}
    /// End function for constructing iterators.
    inline const DataIndex& end() const
    {
      return m_end;
    };
    /// Checks to see whether DBL is coarsenable. Requires ProblemDomain to be coarsenable, and m_boxSize%a_refrat == Point::Zeros.
    inline bool coarsenable(const Point& a_refrat) const;

    /// Returns DBL coarsened by a_refrat. Checks by PR_assert whether this is possible.
    inline DisjointBoxLayout coarsen(const Point& a_refrat) const;

    /// Returns DBL refined by a_refrat. 
    inline DisjointBoxLayout refine(const Point& a_refrat) const;

#ifdef PROTO_CUDA
      //Returns the single stream
      static protoStream_t getCurrentStream()
      {
          static protoStream_t currStream;
          static bool init=false;
          if(!init)
          {
              protoStreamCreate(&currStream);
              init=true;
          }
          return currStream;
      }
#endif
  };
  /// Stream output for DBL.
  inline std::ostream& operator<< (std::ostream& os, const DisjointBoxLayout& a_dbl)
  {
    os << "DisjointBoxLayout: \n" ;
    os << a_dbl.problemDomain() << "\n";
    auto internals = a_dbl.dblInternals();
    os << "coarsened Domain = " << internals->coarsenedDomain << "\n" ;
    os << "Points and ranks: \n";
    for (int k = 0; k < internals->allBoxes.size();k++)
      {
        os << internals->allBoxes[k].first 
           << " , " << internals->allBoxes[k].second << "\n";
      }
    os << "StartProc vector: \n";
    for (int k = 0; k < internals->startProc.size() ; k++)
      {
        os << "Starting patch number for processor " << k 
           << " = " << internals->startProc[k] << "\n";
      }
    return os;
  };
#include "implem/Proto_DisjointBoxLayoutImplem.H"
} //end namespace Proto

#endif
