#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif
#ifndef Proto_LEVELITERATORS_H_
#define Proto_LEVELITERATORS_H_
#include <vector>
#include "Proto.H"
#include "Proto_DataIndex.H"
#include "Proto_DisjointBoxLayout.H"

namespace Proto
{
  class DataIterator
  {
  protected:
    DisjointBoxLayout m_dbl;
    DataIndex m_current;
    DataIndex m_end;
  public:
    /// Constructor. default-constructs m_current.
    // DataIterator();
    inline DataIterator(const DisjointBoxLayout& a_dbl);
    //~DataIterator();
    /// Standard iterator member functions begin, end, operator++. 
    inline DataIterator& begin();
    inline DataIndex& end(){return m_end;};
    inline DataIterator& operator++();
    /// Dereference returns a DataIndex.
    inline const DataIndex& operator*() const;
    /// Returns the Point corresponding to the Box for that index. Undefined if invalid.
    inline Point point() const;
    /// Returns the Box corresponding to that index. Undefined if invalid.
    inline Box box() const; 
  };
  class NeighborIterator
  {
  protected:
    DisjointBoxLayout m_dbl;
    Box m_bx;
    BoxIterator m_bxit;
    Point m_shift;
    DataIndex m_current;
    DataIndex m_end;
    ProblemDomain m_dProbC;
  public:
    /// Constructor takes DBL, and Box that we are intersecting with the DBL.
    //NeighborIterator();
    inline NeighborIterator(const DisjointBoxLayout& a_dbl, const Box& a_bx);
    //~NeighborIterator();

    /// Iterator functions begin, end, ++.
    inline NeighborIterator& begin();
    inline const DataIndex& end() const{return m_end;};
    inline NeighborIterator& operator++();
    
    /// Dereference operator returns the DataIndex corresponding to the current iterate.
    inline const DataIndex& operator*() const{return m_current;};
    /// The Point corresponding to the Box for that index. Undefined if invalid.
    inline Point point() const;
    /// The region to be copied from corresponding to that index. Undefined if invalid.
    inline Box srcBox() const;
    /// The region to be copied from corresponding to that index. Undefined if invalid.
    inline Box destBox() const;
    /// The amount by which the srcBox is shifted to obtain destBox a_bx.
    inline Point shift();
    // Processor ID for current iterate.
    inline int procid(){return m_dbl.procid(m_current);}; 
  };
#include "implem/Proto_LevelIteratorsImplem.H"
} // end Proto namespace.
#endif
     
