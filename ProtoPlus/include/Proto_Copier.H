#ifndef __Proto_Copier__
#define __Proto_Copier__

#include "Proto.H"
#include "Proto_DisjointBoxLayout.H"
//  can start using a pool allocation scheme and placement new

namespace proto
{
///
class MotionItem
{
public:
  unsigned int fromIndex, toIndex;
  Box fromRegion;
  Box toRegion;
  unsigned int procID;
  MotionItem(const DataIndex& a_from,
             const DataIndex& a_to,
             const Box&       a_fromRegion,
             const Box&       a_toRegion)
    :fromIndex(a_from), toIndex(a_to), fromRegion(a_fromRegion),
     toRegion(a_toRegion), procID(-1)
  {
  }

  inline bool operator==(const MotionItem& rhs)
  { return fromIndex == rhs.fromIndex && toIndex == rhs.toIndex &&
      fromRegion == rhs.fromRegion && toRegion == rhs.toRegion &&
      procID == rhs.procID;}
};

class CopierBuffer
{
public:

  ///null constructor, copy constructor and operator= can be compiler defined.
  CopierBuffer():m_ncomps(0), m_sendbuffer(NULL), m_sendcapacity(0),
                 m_recbuffer(NULL), m_reccapacity(0)
  {}

  ///
  virtual ~CopierBuffer();

  void clear();

  bool isDefined(int ncomps) const
  { return ncomps == m_ncomps;}

  mutable int m_ncomps;

  mutable void*  m_sendbuffer; // pointer member OK here,
                               // since LevelData<T> has no copy
  mutable size_t m_sendcapacity;
  mutable void*  m_recbuffer;  // pointer member OK here,
                               // since LevelData<T> has no copy
  mutable size_t m_reccapacity;

  struct bufEntry
  {
    void* bufPtr; // pointer into contiguous m_buffer
    size_t size;
    const MotionItem* item;
    unsigned int procID;
  };


  mutable std::vector<bufEntry> m_fromMe;
  mutable std::vector<bufEntry> m_toMe;
  mutable std::vector<std::vector<bufEntry> > m_toMeUnpack;

#ifdef PR_MPI
  mutable Vector<MPI_Request>  m_sendRequests,  m_receiveRequests;
  mutable Vector<MPI_Status>   m_receiveStatus, m_sendStatus;
#endif
  mutable int numSends, numReceives;

protected:

};

/// A strange but true thing to make copying from one boxlayoutdata to another fast
/**
  class to handle the organization of data from a LevelData to a BoxLayoutData
  defined over the same index space.  The idea behind this object
  is twofold:
  A)  there is a significant amount of information that can be
  computed once and reused from one copy to the next for
  a LevelData and BoxLayoutData that have unchanging BoxLayouts
  underneath them.  In these cases, LevelData::copyTo methods
  has an optional argument.  That way, in cases where the operation
  is only performed a single time, we do not waste our time building
  optimizing data, etc.
  B)  when we interface with KeLP, this class maps quite well
  with the existing KeLP API (class MotionPlan).

  pains have been taken to not have this class be templated.  Thus
  it does not actually perform the copying itself (which would require
  knowledge of the template class type).  It provides an iterator
  interface to MotionItems that specify the required work.  LevelData
  provides the copyTo methods.
*/

class Copier
{
public:

  ///null constructor, copy constructor and operator= can be compiler defined.
  Copier():m_isDefined(false)
  {}

  Copier(const Copier& a_rhs);

  ///
  Copier(const DisjointBoxLayout& a_level,
         const BoxLayout& a_dest,
         bool a_exchange = false,
         IntVect a_shift = IntVect::Zero);

  /// this constructor contains support for periodic BC's
  Copier(const DisjointBoxLayout& a_level,
         const BoxLayout& a_dest,
         const ProblemDomain& a_domain,
         bool a_exchange  = false,
         IntVect a_shift = IntVect::Zero);

  /// Copier to copy data into the valid and invalid regions of a_dest from a_level
  Copier(const DisjointBoxLayout& a_level,
         const BoxLayout& a_dest,
         const IntVect& a_destGhost,
         bool a_exchange = false,
         IntVect a_shift = IntVect::Zero);

  /// Copier to copy data into the valid and invalid regions of a_dest from a_level with support for periodic BCs
  Copier(const DisjointBoxLayout& a_level,
         const BoxLayout& a_dest,
         const ProblemDomain& a_domain,
         const IntVect& a_destGhost,
         bool  a_exchange = false,
         IntVect a_shift = IntVect::Zero);

  ///
  virtual ~Copier();

  ///
  virtual void define(const DisjointBoxLayout& a_level,
                      const BoxLayout& a_dest,
                      bool  a_exchange = false,
                      IntVect a_shift = IntVect::Zero);

  /// contains support for periodic BCs
  virtual void define(const DisjointBoxLayout& a_level,
                      const BoxLayout& a_dest,
                      const ProblemDomain& a_domain,
                      bool  a_exchange = false,
                      IntVect a_shift = IntVect::Zero);

  ///
  virtual void define(const DisjointBoxLayout& a_level,
                      const BoxLayout& a_dest,
                      const IntVect& a_destGhost,
                      bool  a_exchange = false,
                      IntVect a_shift = IntVect::Zero);

  /// contains support for periodic BCs
  virtual void define(const BoxLayout& a_level,
                      const BoxLayout& a_dest,
                      const ProblemDomain& a_domain,
                      const IntVect& a_destGhost,
                      bool  a_exchange = false,
                      IntVect a_shift = IntVect::Zero);

  ///  alternative version of define for exchange copiers that uses new optimized neighborIterator
  virtual void exchangeDefine(const DisjointBoxLayout& a_grids,
                              const IntVect& a_ghost,
	                      bool a_includeSelf=false);
  
  void defineFixedBoxSize(const DisjointBoxLayout& a_src,
                          const LMap&  a_lmap,
                          const IntVect&  a_ghost,
                          const ProblemDomain& a_domain,
                          bool  a_includeSelf=false,
                          bool  a_reverse=false);

  void defineFixedSizeNodesCollect(const DisjointBoxLayout& a_layout,
                                   const LMap&  a_lmap,
                                   const IntVect& a_ghostSrc,
                                   const IntVect& a_ghostDst,
                                   const ProblemDomain& a_domain);
  
  /// performs the computation necessary for moving from ghost+valid cells to all valid cells in a_dest
  /** instead of copying from the valid of a_src to the valid+ghost of a_dest,
      this routine defines a Copier that moves data from the ghost+valid region
      of a_src to  a_dest.  This is used by LevelFluxRegister.
  */
  void ghostDefine(const DisjointBoxLayout& a_src,
                   const DisjointBoxLayout& a_dest,
                   const ProblemDomain& a_domain,
                   const IntVect& a_srcGhost);

  /// performs the computation necessary for moving from ghost+valid cells to all valid+ghost cells in a_dest
  /** instead of copying from the valid of a_src to the valid+ghost of a_dest,
      this routine defines a Copier that moves data from the ghost+valid region
      of a_src to valid+ghost a_dest.  This is used by BlockBaseRegister.  Note that this does not handle 
      periodic matchings here.  That is handled by different code in the MappedMultiBlock libraries.
      Sending in a null-constructed ProblemDomain just results in a computation that ignores domain boundaries.
  */
  void ghostDefine(const DisjointBoxLayout& a_src,
                   const DisjointBoxLayout& a_dest,
                   const ProblemDomain& a_domain,
                   const IntVect& a_srcGhost,
                   const IntVect& a_destGhost);

  ///
  virtual void clear();

  ///
  /**
     The assignment operator.

  */
  Copier& operator= (const Copier& b);

  ///
  /**
     reverse the communication of this Copier.
  */
  void reverse();

  ///
  /**
    Under the right algorithmic circumstances, a user might know that they no not
    need corner cells filled in exchanges, what they really need are just sufficient
    ghost cells in each CH_SPACEDIM direction for each valid cell.  This operation
    eliminates copy operations that do not contribute to such cells.  Under reasonable
    3D configurations this can eliminate more than half the copy operations.
  */
  void trimEdges(const DisjointBoxLayout& a_exchangedLayout, const IntVect& ghost);

  ///
  /*
    OK, if you *really* know what you are doing (!)  you can coarsen a Copier in place of
    building a new Copier on a coarsened DisjointBoxLayout.  The trick is when you do have
    ghost cells and working out on paper for yourself that the coarsening is valid. One
    common configuration is an exchange Copier with *one* ghost cell.  A particular
    case where coarsen will do the *wrong* thing is an exchange Copier with more than one
    ghost cell.
  **/
  void coarsen(int a_refRatio);

  ///
  /*  check equality of two Copiers.  Full check of every MotionItem
   */
  bool operator==(const Copier& rhs) const;
  
  int print() const;

  bool bufferAllocated() const;
  void setBufferAllocated(bool arg) const;

  int numLocalCellsToCopy() const;
  int numFromCellsToCopy() const;
  int numToCellsToCopy() const;

  bool isDefined() const
  { return m_isDefined;}

  CopierBuffer  m_buffers;

  std::vector<IndexTM<int,2> >  m_range;
protected:



  Vector<MotionItem*>   m_localMotionPlan;
  Vector<MotionItem*> m_fromMotionPlan;
  Vector<MotionItem*> m_toMotionPlan;

  friend void dumpmemoryatexit();
  static Pool s_motionItemPool;
  mutable bool buffersAllocated;

  // keep a refcounted reference around for debugging purposes, we can
  // decide afterwards if we want to eliminate it.
  DisjointBoxLayout m_originPlan;
  BoxLayout  m_dest;

  bool m_isDefined;

  void trimMotion(const DisjointBoxLayout& a_exchangedLayout, const IntVect& a_ghost,
                  const Vector<MotionItem*>& a_oldItems, Vector<MotionItem*>& a_newItems);

  void sort();

  // sneaky end-around to problem of getting physDomains in derived classes
  const ProblemDomain& getPhysDomain(const DisjointBoxLayout& a_level) const;
};



}
#endif
