#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

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
    DataIterator(const DisjointBoxLayout& a_dbl);
    ~DataIterator();
    /// Standard iterator member functions begin, end, operator++. 
    inline DataIterator& begin();
    inline DataIndex& end(){return m_end;};
    inline DataIterator& operator++();
    /// Dereference returns a DataIndex.
    inline const DataIndex& operator*() const {return m_current;};
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
  public:
    /// Constructor takes DBL, and Box that we are intersecting with the DBL.
    NeighborIterator(const DisjointBoxLayout& a_dbl, const Box& a_bx);
    ~NeighborIterator();

    /// Iterator functions begin, end, ++.
    inline NeighborIterator& begin();
    inline DataIndex& end(){return m_end;};
    inline NeighborIterator& operator++();
    
    /// Dereference operator returns the DataIndex corresponding to the current iterate.
    inline DataIndex operator*(){return m_current;};
    /// The Point corresponding to the Box for that index. Undefined if invalid.
    inline Point point();
    /// The Box corresponding to that index. Undefined if invalid.
    inline Box box();
    /// The amount by which the current Box in the DBL is shifted to intersect a_bx.
    inline Point shift();
    // Processor ID for current iterate.
    inline int procID(){return m_dbl.procID(m_current);}; 
  };
#if 0
  /// Class that manages iteration over the flux register regions for a fine Box.
  /** 
      Used to to construct Copier for flux registers (Copier::FluxRegisterDefine).
      There is potentially one Coarse Box for each (LoHiSide,direction), that exists 
      only if the region defined by the register is not covered by another Box in the 
      fine DBL.
   */
  class FineFRIterator
  {
  public:
    /// Constructor takes DataIndex associated with the fine Box, and the refinement ratio.
   FineFRIterator(const DataIndex& a_fineIndex,
                  const Point& a_refRatio);
    ~FineFRIterator();

    /// Iterator functions begin, end, ++.
    inline void begin();
    inline DataIndex end();
    inline void operator++();
    /// Dereference operator returns the DataIndex corresponding to this object.
    inline DataIndex operator*{return m_fineIndex;};

    /// The low corner of the coarse Box corresponding to the current register.
    inline Point operator();

    /// The coarse Box corresponding to the current register.
    inline Box box();

    /// The shift applied to the coarse Box to get it into the valid region. 
    inline Point srcShift();
  protected:
    DataIndex m_fineIndex;
    pair<int,int> m_currentRegister;
    //first element is 0 or 1 (LoHiSide), second element is dir.
  };
  /// Class that manages iteration over the flux register regions for a coarse Box.
  class CoarseFRIterator
  {
  public:
    /// Coanstructor. 
    CoarseFRIterator(const& a_DisjointBoxLayout m_fineDBL,
                     const DataIndex& a_coarseDit,
                     const Point& a_refRatio);
     ~CoarseFRIterator();
    void begin();
    bool end();
    void operator++();
    DataIndex operator();
    Box box();
    Point srcShift();
    DataIndex fromIndex();
  protected:
    DataIndex m_coarseDit;
    DisjointBoxlayout m_fine;
    pair<LoHiSide,int> m_currentRegisterSet;
    BoxIterator m_bxit;
    // For each LoHiSide, direction, we iterate over a box corresponding to the
    // fine boxes coarsened to single Points inside the coarse box.
  }
} // end Proto namespace.
#include "Proto_LevelIteratorsImplem.H"
#endif
