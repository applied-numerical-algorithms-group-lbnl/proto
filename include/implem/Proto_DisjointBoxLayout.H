#ifndef _PROTO_DisjointBoxLayout_H___
#define _PROTO_DisjointBoxLayout_H___

#include "Proto_Point.H"
#include "Proto_Box.H"
#include "Proto_Stack.H"
#include <cstdlib> //for size_t
#include <iostream>
#include <stack>
#include <memory>

using std::shared_ptr;
using std::array;
using std::vector;
using std::stack;
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
#ifdef PROTO_CUDA


    struct cudaData
    {
      int s_numStreams;
      std::vector<StackAlloc>   s_allocator;
      std::vector<cudaStream_t>::size_type s_currentStream;
      std::vector<cudaStream_t> s_streams;
    };

    
    static shared_ptr<cudaData> getCudaData() 
    {
      static shared_ptr<cudaData>  s_cudaData;
      static bool alloced = false;
      int numstreams = 8;
      if(!alloced)
      {
        s_cudaData = shared_ptr<cudaData>(new cudaData());
        
        s_cudaData->s_numStreams    = numstreams;
        s_cudaData->s_currentStream = 0;
        alloced = true;
      }
      return s_cudaData;
    }

    static void setNumStreams(int a_numstreams) 
    {
      static shared_ptr<cudaData> cudat = getCudaData();
      cudat->s_numStreams    = a_numstreams;
      initStreams();
      
    }
    static void activeBox(unsigned int a_id)
    {
      PR_TIME("activeBox");
      static shared_ptr<cudaData> cudat = getCudaData();
      if((a_id - cudat->s_numStreams)%cudat->s_numStreams != cudat->s_currentStream)
        {
          cudat->s_currentStream = (++cudat->s_currentStream)%(cudat->s_numStreams);
          clearStream(cudat->s_currentStream);
        }
    }
    static cudaStream_t getCurrentStream()
    {
      static shared_ptr<cudaData> cudat = getCudaData();
      //this can get called outside context of DBL so make certain stuff got initialized
      initStreams(); //null op if alread initialized
      return cudat->s_streams[cudat->s_currentStream];
    }

    static StackAlloc& getStackAllocator()
    {
      static shared_ptr<cudaData> cudat = getCudaData();
      //this can get called outside context of DBL so make certain stuff got initialized
      initStreams(); //null op if alread initialized
      return cudat->s_allocator[cudat->s_currentStream];
    }

    static  void clearStream(unsigned int streamID)
    {
      PR_TIME("clearStream");
      static shared_ptr<cudaData> cudat = getCudaData();
      cudaStreamSynchronize(cudat->s_streams[streamID]);
      cudat->s_allocator[cudat->s_currentStream].reset();
    }
      
    static bool initStreams()
    {
      static bool initialized = false;
      static shared_ptr<cudaData> cudat = getCudaData();
      if(!initialized || (cudat->s_streams.size() != cudat->s_numStreams))
      {
        initialized = true; 
        if(cudat->s_streams.size() != cudat->s_numStreams)
        {
          for(std::vector<cudaStream_t>::size_type i=0; i< cudat->s_streams.size(); i++)
          {
            cudaStreamDestroy(cudat->s_streams[i]);
          }
          cudat->s_streams.resize(cudat->s_numStreams);
          cudat->s_allocator.resize(cudat->s_numStreams);
          for(std::vector<cudaStream_t>::size_type i=0; i< cudat->s_streams.size(); i++)
          {
            cudaStreamCreate(&cudat->s_streams[i]);
          }
        }
      }
      return true;
    }

#else
    struct hostData
    {
      StackAlloc   s_allocator;
      unsigned int s_current;
    };
    static shared_ptr<hostData> getHostData() 
    {
      static shared_ptr<hostData>  s_hostData;
      static bool alloced = false;
      if(!alloced)
      {
        s_hostData = shared_ptr<hostData>(new hostData());
        s_hostData->s_current = 0;
        alloced = true;
      }
      return s_hostData;
    }

    static StackAlloc& getStackAllocator()
    {
      static shared_ptr<hostData> hostdat = getHostData();
      return hostdat->s_allocator;
    }
    static  void clearMemory()
    {
      static shared_ptr<hostData> hostdat = getHostData();
      hostdat->s_allocator.reset();
    }
    static void activeBox(unsigned int a_id)
    {
      static shared_ptr<hostData> hostdat = getHostData();
      if(a_id != hostdat->s_current)
      {
        clearMemory();
        hostdat->s_current = a_id;
      }
    }

#endif
  };
} //end namespace Proto
#include "Proto_DisjointBoxLayoutI.H"
#endif
