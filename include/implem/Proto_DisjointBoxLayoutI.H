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
  {
#ifdef PROTO_CUDA
    initStreams();
#endif
  }

  ///
  DisjointBoxLayout::
  DisjointBoxLayout(const Box             & a_domain, 
                    const unsigned int    & a_maxgrid,
                    const array<bool,DIM> & a_periodic)
  {
#ifdef PROTO_CUDA
    initStreams();
#endif
    define(a_domain, a_maxgrid, a_periodic);
  }

  ///
  DisjointBoxLayout::
  DisjointBoxLayout(const DisjointBoxLayout& a_input)
  {
#ifdef PROTO_CUDA
    initStreams();
#endif
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



/**********/
  inline Point getPeriodicShift(Point a_neigh, Box a_coarDom)
  {
    //this routine is to find where we have to shift to get back into the domain
    PR_assert(!a_coarDom.contains(a_neigh));
    
    Point dcoarLo = a_coarDom.low();
    Point dcoarHi = a_coarDom.high();
    //once we are outside the domain for a direction,
    //do the cases to figure out the periodic shift
    Point perShift = Point::Zeroes();
    for(int idir = 0; idir < DIM; idir++)
    {
      
      int  signedShiftI = 0; //initialization here is important
      if(     a_neigh[idir] < dcoarLo[idir])
      {
        signedShiftI =  a_coarDom.size(idir);
      }
      else if(a_neigh[idir] > dcoarHi[idir])
      {
        signedShiftI = -a_coarDom.size(idir);
      }
    
      perShift += signedShiftI*Point::Basis(idir);
    }

    return perShift;
  }
  ///useful for copyto and exchange
  vector<NeighborElement> 
  DisjointBoxLayout::
  getNeighbors(unsigned int a_index)
  {
    vector<NeighborElement> retval;
    const Box &  coarDom = m_internals->m_coarsenedDom;
    Point        coarPt  = m_internals->m_coarsenedDom[a_index];
    Box  coarBox(coarPt, coarPt);
    Box neighBox = coarBox.grow(1);

    for(unsigned int ipt = 0; ipt < neighBox.size(); ipt++)
    {
      Point neighPtCoar = neighBox[ipt];
      //avoid self copy
      if(coarPt != neighPtCoar)
      {
        NeighborElement neighElem;        
        if(coarDom.contains(neighPtCoar))
        {
          neighElem.m_periodicShift  = Point::Zeros();
          neighElem.m_boxID = coarDom.index(neighPtCoar);
          retval.push_back(neighElem);
        }
        else
        {
          Point perShiftCoar = getPeriodicShift(neighPtCoar, coarDom);
          neighPtCoar += perShiftCoar;
          if(coarDom.contains(neighPtCoar))  //to avoid weird more ghost than domain bugs
          {
            neighElem.m_periodicShift  = perShiftCoar*(m_internals->m_maxgrid);
            neighElem.m_boxID = coarDom.index(neighPtCoar);
            retval.push_back(neighElem);
          }

        }
      }
    }

    return retval;
  }

  //this does the corners
} //end namespace Proto
#endif
