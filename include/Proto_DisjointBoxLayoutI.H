#ifndef _PROTO_DisjointBoxLayoutImplem_H_
#define _PROTO_DisjointBoxLayoutImplem_H_

#include "Proto_DisjointBoxLayout.H"
#include <cstdlib> //for size_t

using std::shared_ptr;
using std::array;
namespace Proto 
{

  DisjointBoxLayout::
  DisjointBoxLayout()
  {;}

  ///
  DisjointBoxLayout::
  DisjointBoxLayout(const Box             & a_domain, 
                    const unsigned int    & a_maxgrid,
                    const array<bool,DIM> & a_periodic)
  {
    define(a_domain, a_maxgrid, a_periodic);
  }

  ///
  DisjointBoxLayout::
  DisjointBoxLayout(const DisjointBoxLayout& a_input)
  {
    if(&a_input != this)
    {
      m_internals = a_input.m_internals;
    }
  }

  ///
  DisjointBoxLayout& 
  DisjointBoxLayout::operator=(const DisjointBoxLayout& a_input)
  {
    if(&a_input != this)
    {
      m_internals = a_input.m_internals;
    }
    return *this;
  }

  ///
  bool 
  DisjointBoxLayout::
  operator==(const DisjointBoxLayout& a_input) const
  {
    return (m_internals == a_input.m_internals);
  }

  ///
  void 
  DisjointBoxLayout::
  define(const Box             & a_domain, 
         const unsigned int    & a_maxgrid,    
         const array<bool,DIM> & a_periodic)
  {
    PROTO_ASSERT(a_domain.coarsenable(a_maxgrid), "invalid dbl combo");

    m_internals = shared_ptr<localData>(new localData());

    m_internals->m_inputDomain  = a_domain;
    m_internals->m_coarsenedDom = a_domain.coarsen(a_maxgrid);
    m_internals->m_maxgrid      = a_maxgrid;
    m_internals->m_periodic     = a_periodic;
  }

  ///
  Box 
  DisjointBoxLayout::
  operator[](unsigned int a_index) const
  {
    PROTO_ASSERT(m_internals,"trying to access undefined dbl boxes");
    Point coarpt = m_internals->m_coarsenedDom[a_index];
    Box coarBox(coarpt, coarpt);
    Box retval = coarBox.refine(m_internals->m_maxgrid);
    return retval;
  }


  ///number of boxes in grid (over all procs)
  unsigned int 
  DisjointBoxLayout::
  size() const
  {
    return m_internals->m_coarsenedDom.size();
  }
  
  ///useful for copyto and exchange
  vector<NeighborElement> 
  DisjointBoxLayout::
  getNeighbors(unsigned int a_index)
  {
    vector<NeighborElement> retval;
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
        neighElem.m_periodicShift  = perShiftCoar*(m_internals->m_maxgrid);
        neighElem.m_boxID = coardom.index(neighPtCoar);
      }
      retval.push_back(neighElem);
    }

    return retval;
  }

}; //end namespace Proto
#endif
