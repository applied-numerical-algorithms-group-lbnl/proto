#ifndef __Proto_DisjointBoxLayout__
#define __Proto_DisjointBoxLayout__
#include <cstdio>
#include <cstring>
#include <cassert>
#include <cmath>

#include <vector>
#include <memory>

#include <iostream>
#include <fstream>
#include <sstream>


#include "Proto.H"
namespace proto
{


///boxes that cover a domain box
  class DisjointBoxLayout
  {
  private:
    struct localData
    {
      Box                    m_coarsenedDom;
      vector< unsigned int>  m_localBoxes;
      unsigned int           m_maxgrid;
      vector< unsigned int>  m_procs;
    };

    //this is to make this a ref-counted object
    shared_ptr<localData> m_internals;


  public:

    DisjointBoxLayout()
    {;}

    ///
    DisjointBoxLayout(const Box& a_domain, const unsigned int& a_maxgrid)
    {
      define(a_domain, a_maxgrid);
    }


    ///
    DisjointBoxLayout(const DisjointBoxLayout& a_input)
    {
      if(&a_input != this)
      {
        m_internals = a_input.m_internals;
      }
    }

    ///
    DisjointBoxLayout& operator=(const DisjointBoxLayout& a_input)
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
    inline void define(const Box& a_domain, const unsigned int& a_maxgrid)
    {
      PROTO_ASSERT(a_domain.coarsenable(a_maxgrid), "invalid dbl combo");

      m_internals = shared_ptr<localData>(new localData());

      m_internals->m_coarsenedDom = a_domain.coarsen(a_maxgrid);
      m_internals->m_maxgrid = a_maxgrid;
      m_internals->m_procmap.define(m_coarsenedDom);

      //should probably do some sort of nearest neighbor walk
      unsigned int numboxes = m_internals->m_coarsenedDom.size();
      m_internals->m_procs.resize(numboxes);
      unsigned int boxesperproc = numboxes/(CH_numProc());
      for(unsigned int ibox = 0; ibox < numboxes; ibox++)
      {
        unsigned int boxproc = ibox/boxesperproc;
        m_internals->m_procs[ibox] = boxproc;
        unsigned int procid = CH_procID();
        if(boxproc == procid)
        {
          m_internals->m_localBoxes.push_back(ibox);
        }
      }

    }

    ///
    inline unsigned  int procID(unsigned int a_index) const
    {
      PROTO_ASSERT(m_internals,"trying to access undefined dbl procids");
      return m_internals->m_procs[a_index];
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


    ///number of boxes in grid (over all procs)
    inline unsigned int size() const
    {
      return m_internals->m_coarsenedDom.size();
    }

    ///boxes in grid whose data lives on the current proc
    inline const vector<unsigned int>& localBoxes() const
    {
      PROTO_ASSERT(m_internals,"trying to access undefined dbl local boxes");
      return m_internals->m_localBoxes;
    }

  };
}

#endif
