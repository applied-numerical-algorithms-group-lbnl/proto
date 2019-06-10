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

    inline void dump() const
    {
      using std::cout;
      using std::endl;
      cout << "dbl domain   = " << m_internals->m_inputDomain  << endl;
      cout << "dbl coardom  = " << m_internals->m_coarsenedDom << endl;
      cout << "dbl maxgrid  = " << m_internals->m_maxgrid      << endl;
      cout << "dbl boxes: "   << endl;
      for(int ibox = 0; ibox  < size(); ibox++)
      {
        cout << ibox << " " << (*this)[ibox]      << endl;
      }
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

      PROTO_ASSERT(coarsenable(a_refrat), "refrat too big");

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
    static int s_numStreams=8;
#ifdef PROTO_CUDA
    static std::vector<cudaStream_t>::size_type s_currentStream=0;
    static std::vector<cudaStream_t> s_streams;
    static std::vector<std::stack<std::shared_ptr<void> > > s_streamStack;

    static void activeBox(unsigned int a_id)
    {
      if(a_id%s_numStreams != s_currentStream)
        {
          clearStream(s_currentStream);
          s_currentStream ++;
          s_currentStream = s_currentStream%s_numStreams;
        }
    }
    static  void clearStream(unsigned int streamID)
    {
      cudaStreamSynchronize(s_streams[streamID]);
      std::stack<std::shared_ptr<void>>& st = s_streamStack[streamID];
      while(!st.empty())
        {
          st.pop();
        }
    }
      
    static void pushStack(std::shared_ptr<void>& a_data)
    {
      std::stack<std::shared_ptr<void>>& st = s_streamStack[s_currentStream];
      st.push(a_data);
    }
    static bool initStreams()
    {
      if(s_streams.size() != s_numStreams)
        {
          for(s_stream::size_type i=0; i<s_stream.size(); i++)
            {
              cudaStreamDestroy(s_streams[i]);
            }
          s_streams.resize(s_numStreams);
          for(s_streams::size_type i=0; i<s_streams.size(); i++)
            {
              cudaStreamCreate(&s_streams[i]);
            }
        }
      return true;
    }
  

    }
#else
    static bool initStreams() {return true;}
#endif
  };
} //end namespace Proto
#include "Proto_DisjointBoxLayoutI.H"
#endif
