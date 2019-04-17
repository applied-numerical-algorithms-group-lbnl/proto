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
      Box                   m_inputDomain;
      Box                   m_coarsenedDom;
      unsigned int          m_maxgrid;
      array<bool, DIM>      m_periodic;
    };

    //this is to make this a ref-counted object
    shared_ptr<localData> m_internals;


  public:

    ///
    inline DisjointBoxLayout();

    ///
    inline DisjointBoxLayout(const Box             & a_domain, 
                             const unsigned int    & a_maxgrid,
                             const array<bool,DIM> & a_periodic);


    ///
    inline DisjointBoxLayout(const DisjointBoxLayout& a_input);

    ///
    inline DisjointBoxLayout& operator=(const DisjointBoxLayout& a_input);


    ///
    inline bool operator==(const DisjointBoxLayout& a_input) const;

    ///
    inline void define(const Box             & a_domain, 
                       const unsigned int    & a_maxgrid,
                       const array<bool,DIM> & a_periodic);

    inline Box getDomain() const 
    {
      return m_internals->m_inputDomain;
    }

    ///
    inline Box operator[](unsigned int a_index) const;


    ///number of boxes in grid 
    inline unsigned int size() const;

    ///useful for copyto and exchange
    inline vector<NeighborElement> getNeighbors(unsigned int a_index);

    inline bool coarsenable(int a_refrat) const
    {
      bool retval = m_internals->m_inputDomain.coarsenable(a_refrat);
      if((m_internals->m_maxgrid)%a_refrat != 0) retval = false;
      return retval;
    }

    ///
    inline void coarsen(int a_refrat) 
    {

      PR_assert(coarsenable(a_refrat));

      //redefine with smaller maxgrid
      unsigned int newmaxgrid = m_internals->m_maxgrid/a_refrat;
      Box newdomain = m_internals->m_inputDomain.coarsen(a_refrat);
      define(newdomain, newmaxgrid, m_internals->m_periodic);
    }

    ///
    inline void refine(int a_refrat) 
    {
      unsigned int newmaxgrid = m_internals->m_maxgrid*a_refrat;
      //redefine with bigger maxgrid
      Box newdomain = m_internals->m_inputDomain.refine(a_refrat);
      define(newdomain, newmaxgrid, m_internals->m_periodic);
    }

  };
}; //end namespace Proto
#include "Proto_DisjointBoxLayoutI.H"
#endif
