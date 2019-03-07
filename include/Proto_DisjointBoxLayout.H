#ifndef _PROTO_DisjointBoxLayout_H_
#define _PROTO_DisjointBoxLayout_H_

#include "Proto_Point.H"
#include "Proto_Box.H"
#include <cstdlib> //for size_t
#include <iostream>

using std::shared_ptr;
using std::array;
using std::vector;
namespace Proto 
{

  ///
  struct NeighborElement
  {
    ///id of box which is a neighbor
    unsigned int  m_boxID;

    //needed for periodic (shift neighborbox by -m_periodicshift to get it into the same index space)
    Point m_periodicShift;
    
  };

  ///
  /**
   Simplified version Chombo's DBL that has no concept of MPI processes.   So all boxes are local.
   Also, this is not an AMR thing.   A DBL here covers a given box with boxes of a given size.
  */
  class DisjointBoxLayout
  {
  private:

    struct localData
    {
      Box                   m_coarsenedDom;
      unsigned int          m_maxgrid;
      array<bool, DIM>      m_periodic;
    };

    //this is to make this a ref-counted object
    shared_ptr<localData> m_internals;


  public:

    inline DisjointBoxLayout()
    {;}

    ///
    inline DisjointBoxLayout(const Box             & a_domain, 
                             const unsigned int    & a_maxgrid,
                             const array<bool,DIM> & a_periodic)
    {
      define(a_domain, a_maxgrid, a_periodic);
    }



    ///
    inline DisjointBoxLayout(const DisjointBoxLayout& a_input)
    {
      if(&a_input != this)
      {
        m_internals = a_input.m_internals;
      }
    }


    ///
    inline DisjointBoxLayout& operator=(const DisjointBoxLayout& a_input)
    {
      if(&a_input != this)
      {
        m_internals = a_input.m_internals;
      }
      return *this;
    }


    ///
    inline bool operator==(const DisjointBoxLayout& a_input) const
    {
      return (m_internals == a_input.m_internals);
    }

    ///
    inline void define(const Box             & a_domain, 
                       const unsigned int    & a_maxgrid,
                       const array<bool,DIM> & a_periodic)
    {
      PROTO_ASSERT(a_domain.coarsenable(a_maxgrid), "invalid dbl combo");

      m_internals = shared_ptr<localData>(new localData());

      m_internals->m_coarsenedDom = a_domain.coarsen(a_maxgrid);
      m_internals->m_maxgrid      = a_maxgrid;
      m_internals->m_periodic     = a_periodic;
    }



    ///
    inline Box operator[](unsigned int a_index) const
    {
      PROTO_ASSERT(m_internals,"trying to access undefined dbl boxes");
      Point coarpt = m_internals->m_coarsenedDom[a_index];
      Box coarBox(coarpt, coarpt);
      Box retval = coarBox.refine(m_internals->m_maxgrid);
      return retval;
    }



    ///number of boxes in grid 
    inline unsigned int size() const
    {
      return m_internals->m_coarsenedDom.size();
    }

    ///useful for copyto and exchange
    inline vector<NeighborElement> getNeighbors(unsigned int a_index)
    {
      vector<NeighborElement> retval;
      unsigned int maxgrid = m_internals->m_maxgrid;
      const Box &  coardom = m_internals->m_coarsenedDom;
      Point        coarPt  = m_internals->m_coarsenedDom[a_index];
      Box  coarBox(coarPt, coarPt);
      Box neighBox = coarBox.grow(1);
      for(unsigned int ipt = 0; ipt < neighBox.size(); ipt++)
      {
        NeighborElement neighElem;
        Point neighPtCoar = coardom[ipt];
        if(coardom.contains(neighPtCoar))
        {
          neighElem.m_periodicShift  = Point::Zeros();
          neighElem.m_boxID = coardom.index(neighPtCoar);
        }
        else
        {
          Point perShiftCoar = Point::Zeros();
          Point dcoarLo = coardom.low();
          Point dcoarHi = coardom.high();
          //shift to inside the domain
          for(int idir = 0; idir < DIM; idir++)
          {
            if(neighPtCoar[idir] < dcoarLo[idir])
            {
              perShiftCoar[idir] = coardom.size(idir);
            }
            else if(neighPtCoar[idir] < -dcoarHi[idir])
            {
              perShiftCoar[idir] = -coardom.size(idir);
            }
          }
          neighPtCoar += perShiftCoar;
          neighElem.m_periodicShift  = perShiftCoar*maxgrid;
          neighElem.m_boxID = coardom.index(neighPtCoar);
        }
        retval.push_back(neighElem);
      }

      return retval;
    }

  };
#include "Proto_DisjointBoxLayoutI.H"
} //end namespace Proto
#endif
