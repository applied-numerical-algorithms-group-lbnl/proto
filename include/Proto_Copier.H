#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _PROTO_COPIER_H_
#define _PROTO_COPIER_H_

#include "Proto_DisjointBoxLayout.H"
#include "Proto.H"
#include "Proto_SPMD.H"
#include "Proto_LevelIterators.H"
#include <unordered_map>
#include <cstdint>

namespace Proto
{

  class CopyIterator;

  //  These classes are public because I can't find a nice
  //  way to make a class a friend of all the instantiations
  //  of a template class.  These classes are not part of
  //  the public interface for the Array API.
  //
  //  Later, if MotionItem shows up in the profiler, we
  //  can start using a pool allocation scheme and placement new

  class MotionItem
  {
  public:
    DataIndex fromIndex, toIndex;
    Box fromRegion;
    Box toRegion;
    int procID;
    // this constructor will probably eventually go away
    MotionItem(const DataIndex& a_from,
               const DataIndex& a_to,
               const Box&       a_region);

    MotionItem(const DataIndex& a_from,
               const DataIndex& a_to,
               const Box&       a_fromRegion,
               const Box&       a_toRegion);

    bool operator==(const MotionItem& rhs)
    { return (fromIndex == rhs.fromIndex) && (toIndex == rhs.toIndex) &&
        (fromRegion == rhs.fromRegion) && (toRegion == rhs.toRegion) &&
        (procID == rhs.procID);}
  };

  class CopierBuffer
  {
  public:

    ///null constructor, copy constructor and operator= can be compiler defined.
    CopierBuffer():m_isDefined(false), m_sendbuffer(NULL), m_sendcapacity(0),
                   m_recbuffer(NULL), m_reccapacity(0)
    {}

    ///
    inline virtual ~CopierBuffer();

    inline void clear();

    inline bool isDefined() const
    { return m_isDefined;}

    mutable bool m_isDefined;

    mutable void*  m_sendbuffer; // pointer member OK here,
    // since LevelBoxData has no copy
    mutable size_t m_sendcapacity;
    mutable void*  m_recbuffer;  // pointer member OK here,
    // since LevelBoxData has no copy
    mutable size_t m_reccapacity;

#ifndef DOXYGEN

    struct bufEntry
    {
      void* bufPtr; // pointer into contiguous m_buffer
      size_t size;
      const MotionItem* item;
      unsigned int procID;
      // PC:: need to find out what is going on here.
      bool operator < (const bufEntry& rhs) const
      {
        if (procID == rhs.procID)
          {
            const Box& left = item->toRegion;
            const Box& right= rhs.item->toRegion;
            if (left.low() == right.low())
              {
                return left.high() < (right.high());
              }
            else
              {
                return item->toRegion < rhs.item->toRegion;
              }
          }
        //else
        return procID < rhs.procID;
      }
    };

#endif
    mutable vector<bufEntry> m_fromMe;
    mutable vector<bufEntry> m_toMe;
    mutable vector<vector<bufEntry> > m_toMeUnpack;

//TODO (Aug 24, 2020): These aren't currently used. Remove?
#ifdef PR_MPI
    mutable vector<MPI_Request>  m_sendRequests,  m_receiveRequests;
    mutable vector<MPI_Status>   m_receiveStatus, m_sendStatus;
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
    inline virtual ~Copier();

    ///  alternative version of define for exchange copiers that uses new optimized neighborIterator
    inline virtual void exchangeDefine(const DisjointBoxLayout& a_dbl,
                                const Point& a_ghost);
    /// CopyTo Copier define.
    inline void LBDCopyToDefine(const DisjointBoxLayout& a_src,
                         const DisjointBoxLayout& a_dest);
    ///
    inline virtual void clear();

    ///
    /**
       The assignment operator.

    */
    inline Copier& operator= (const Copier& b);

    ///
    /**
       reverse the communication of this Copier.
    */
    inline void reverse();

    ///
    /*  check equality of two Copiers.  Full check of every MotionItem
     */
    inline bool operator==(const Copier& rhs) const;
  
    inline int print() const;

    inline bool bufferAllocated() const;
    inline void setBufferAllocated(bool arg) const;

    inline int numLocalCellsToCopy() const;
    inline int numFromCellsToCopy() const;
    inline int numToCellsToCopy() const;

    inline bool isDefined() const
    { return m_isDefined;}

    CopierBuffer  m_buffers;


  protected:

    friend class CopyIterator;


    vector<MotionItem*> m_localMotionPlan;
    vector<MotionItem*> m_fromMotionPlan;
    vector<MotionItem*> m_toMotionPlan;

    inline friend void dumpmemoryatexit();
    mutable bool buffersAllocated;

    // keep a refcounted reference around for debugging purposes, we can
    // decide afterwards if we want to eliminate it.
    DisjointBoxLayout m_originPlan;
    DisjointBoxLayout  m_dest;

    bool m_isDefined;

    inline void sort();

  };

  inline std::ostream& operator<<(std::ostream& os, const Copier& copier);

  //===========================================================================
  // end of public interface for Copier.
  //===========================================================================

  inline MotionItem::MotionItem(const DataIndex& a_from,
                                const DataIndex& a_to,
                                const Box&       a_region)
    :fromIndex(a_from), toIndex(a_to), fromRegion(a_region),
     toRegion(a_region), procID(-1)
  {
  }

  inline MotionItem::MotionItem(const DataIndex& a_from,
                                const DataIndex& a_to,
                                const Box&       a_fromRegion,
                                const Box&       a_toRegion)
    :fromIndex(a_from), toIndex(a_to), fromRegion(a_fromRegion),
     toRegion(a_toRegion), procID(-1)
  {
  }
  class CopyIterator
  {
  public:
    enum local_from_to
      {
        LOCAL,
        FROM,
        TO
      };

    inline CopyIterator(const Copier& a_copier, local_from_to);

    inline const MotionItem& operator()() const;

    inline const MotionItem& operator[](size_t index) const;

    inline void operator++();

    inline bool ok() const;

    inline void reset();

    inline size_t size(){ return m_motionplanPtr->size();}

  private:
    const vector<MotionItem*>* m_motionplanPtr;
    unsigned int  m_current;
  };

  //====== inlined functions====================================

  inline CopyIterator::CopyIterator(const Copier& a_copier, local_from_to type)
    :m_current(0)
  {
    switch(type)
      {
      case LOCAL:
        m_motionplanPtr = &(a_copier.m_localMotionPlan);
        break;
      case FROM:
        m_motionplanPtr = &(a_copier.m_fromMotionPlan);
        break;
      case TO:
        m_motionplanPtr = &(a_copier.m_toMotionPlan);
        break;
        //PC:: huh ? default:
        //PC::MayDay::Error("illegal local_from_to option for CopyIterator");
      }
  }

  inline const MotionItem& CopyIterator::operator()() const
  {
    PR_assert(m_current < m_motionplanPtr->size());
    return *(m_motionplanPtr->operator[](m_current));
  }

  inline const MotionItem& CopyIterator::operator[](size_t a_index) const
  {
    PR_assert(a_index < m_motionplanPtr->size());
    return *(m_motionplanPtr->operator[](a_index));
  }

  inline void CopyIterator::operator++()
  {
    ++m_current;
  }

  inline bool CopyIterator::ok() const
  {
    return m_current < m_motionplanPtr->size();
  }

  inline void CopyIterator::reset()
  {
    m_current = 0;
  }
} // end namespace Proto.
#include "implem/Proto_CopierImplem.H"
#endif
