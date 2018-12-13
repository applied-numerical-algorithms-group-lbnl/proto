#include <cstdio>
#include <cstring>
#include <cassert>
#include <cmath>

#include <vector>
#include <memory>

#include <iostream>
#include <fstream>
#include <sstream>

#define DIM 2
#include "../../include/Proto.H"
using std::cout;
using std::endl;
using std::vector;
using std::shared_ptr;
using namespace Proto;
/**/
void
parseCommandLine(unsigned int & a_nx, unsigned int& a_maxgrid, unsigned int & a_numapplies, int argc, char* argv[])
{
  //defaults
  a_nx = 16;
  a_numapplies = 1;
  a_maxgrid = 4;
  cout << "kernel timings of various laplacians" << endl;
  cout << "usage:  " << argv[0] << " -n nx[default:8] -a num_iterations[default:1]" << endl;
  for(int iarg = 0; iarg < argc-1; iarg++)
  {
    if(strcmp(argv[iarg],"-n") == 0)
    {
      a_nx = atoi(argv[iarg+1]);
    }
    else if(strcmp(argv[iarg], "-a") == 0)
    {
      a_numapplies = atoi(argv[iarg+1]);
    }
    else if(strcmp(argv[iarg], "-m") == 0)
    {
      a_maxgrid = atoi(argv[iarg+1]);
    }
  }
}


///proxies for Chombo-style SPMD functions
unsigned int numProc()
{
  return 1;
}

unsigned int procID()
{
  return 0;
}

///boxes that cover a domain box
class DisjointBoxLayout
{
private:
  class localData
  {
  public:
    localData() {;}
    Box                   m_coarsenedDom;
    vector<unsigned int>  m_procs;
    vector<unsigned int>  m_localBoxes;
    unsigned int          m_maxgrid;
  };

  //this is to make this a ref-counted object
  shared_ptr<localData> m_internals;


public:
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
      m_internals = a_input.internals;
    }
  }

  ///
  DisjointBoxLayout& operator=(const DisjointBoxLayout& a_input)
  {
    if(&a_input != this)
    {
      m_internals = a_input.internals;
    }
    return *this;
  }

  ///
  bool operator==(const DisjointBoxLayout& a_input) const
  {
    return (m_internals == a_input.m_internals);
  }

  ///
  void define(const Box& a_domain, const unsigned int& a_maxgrid)
  {
    PROTO_ASSERT(a_domain.coarsenable(a_maxgrid), "invalid dbl combo");

    m_internals = shared_ptr<localData>(new localData());

    m_internals->m_coarsenedDom = a_domain.coarsen(a_maxgrid);
    m_internals->m_maxgrid = a_maxgrid;

    //should probably do some sort of nearest neighbor walk
    unsigned int numboxes = m_internals->m_coarsenedDom.size()
    m_internals->m_procs.resize(numboxes);
    unsigned int boxesperproc = numboxes/numProc();
    m_internals->m_numBoxesThisProc = 0;
    for(unsigned int ibox = 0; ibox > numboxes; ibox++)
    {
      unsigned int boxproc = ibox/boxesperproc;
      m_internals->m_procs[ibox] = boxproc;
      if(boxproc == procID())
      {
        m_internals->m_localBoxes.push_back(ibox);
      }
    }

  }

  ///
  unsigned  int procID(int a_index) const
  {
    PROTO_ASSERT(m_internals,"trying to access undefined dbl procids");
    Point coarpt = m_internals->m_coarsenedDom.index(a_index);
    return m_internals->m_procs[a_index];
  }

  ///
  Box operator[](int a_index) const
  {
    PROTO_ASSERT(m_internals,"trying to access undefined dbl boxes");
    Point coarpt = m_internals->m_coarsenedDom.index(a_index);
    Box retval(coarpt, coarpt);
    retval.refine(m_internals->m_maxgrid);
    return retval;
  }


  ///number of boxes in grid (over all procs)
  unsigned int size() const
  {
    return m_internals->m_coarsenedDom.size();
  }

  ///boxes in grid whose data lives on the current proc
  const vector<unsigned int>& localBoxes() const
  {
    PROTO_ASSERT(m_internals->m_isDefined,"trying to access undefined dbl size");
    return m_internals->m_localBoxes;
  }

};

///data over a disjointboxlayout with ghost cells
template <class T>
class LevelData
{

  ///get to the data on a particular box.  this index is into m_data---the boxes on THIS processor.
  /**
     you can get a vector of these boxes by calling DisjointBoxLayout::localBoxes
   */
  T & operator[](unsigned int a_index)
  {
    
    PROTO_ASSERT(m_isDefined,"trying to access undefined leveldata");
    PROTO_ASSERT(a_index < m_data.size(),"bogus index sent to leveldata");
    
    return (*(m_data[a_index]));
  }

  ///
  LevelData()
  {
    m_isDefined = false;
  }


  ///
  LevelData(const DisjointBoxLayout& a_grids)
  {
    define(a_grids);
  }

  ///
  void define(const DisjointBoxLayout& a_grids, const Point& a_grow)
  {
    m_isDefined = true;
    m_grids = a_grids;
    const vector<unsigned int>& localBoxes = a_grids.localBoxes();
    m_data.resize(localBoxes.size());
    for(unsigned int ibox = 0; ibox < localBoxes.size(); ibox++)
    {
      int boxid   = localBoxes[ibox];
      Box dblbox  = a_grids[boxid];
      Box databox = dblbox.grow(a_grow);
      m_data[ibox] = shared_ptr<T>(new T(databox));
    }
    
  }

  ///return the number of boxes on THIS proc
  unsigned int size() const
  {
    return m_data.size();
  }

  void setVal(const T& a_val)
  {
    for(unsigned int ibox = 0; ibox < m_data.size(); ibox++)
    {
      m_data[ibox].setVal(a_val);
    }
  }
private: 
  //in parallel, this is be the data on this proc
  vector<shared_ptr<T> >        m_data;
  DisjointBoxLayout             m_grids;
  bool                          m_isDefined();
};


///
inline void sync()
{
  #ifdef PROTO_CUDA
    {
      PR_TIME("device sync");
      cudaDeviceSynchronize();
    }
#endif
}
/**/

template <class T> void
applyLaplacians(LevelData< BoxData<T> > & a_phi,
                LevelData< BoxData<T> > & a_lap,
                const DisjointBoxLayout & a_dbl,
                const unsigned int      & a_numapplies)
{

  PR_TIME("applyLaplacians");
#if DIM==2
  Stencil<T> lapSten = Stencil<T>::Laplacian_9();
#else 
  Stencil<T> lapSten = Stencil<T>::Laplacian_27();
#endif

  //remember this is just for timings
  a_phi.setVal(0.);
  a_lap.setVal(0.);
  vector<unsigned int> localBoxes = a_dbl.localBoxes();
  vector<cudaStream_t> streams(localBoxes.size());
  for(unsigned int ibox = 0; ibox < localBoxes.size(); ibox++)
  {
    cudaStreamCreate(&streams[ibox]);
  }
  {
    PR_TIME("applyLaplacianStencilOnLevel");
    for(unsigned int iapp = 0; iapp < a_numapplies; iapp++)
    {
      for(unsigned int ibox = 0; ibox < localBoxes.size(); ibox++)
      {
        Box appBox       = a_dbl[localBoxes[ibox]];
        int streamNumber = ibox;  //perhaps this needs to be thought out better
        lapSten.cudaApplyStream(phi[ibox], lap[ibox], appBox, true, 1.0, streams[ibox]);
      }
    }
    sync();
  }

//  Stencil<T> emptySten;
}
/**/
int main(int argc, char* argv[])
{
  //have to do this to get a time table
  PR_TIMER_SETFILE("proto.time.table");
  unsigned int nx, niter, maxgrid;
  parseCommandLine(nx, maxgrid, niter, argc, argv);

  Point lo = Point::Zeros();
  Point hi = Point::Ones(nx - 1);
  Box domain(lo, hi);
  
  DisjointBoxLayout dbl(domain, maxgrid);
  LevelData<BoxData<double> > phid, lapd;
  LevelData<BoxData<float>  > phif, lapf;

  {
    
    PR_TIME("data definition");

    Point ghostPt = Point::Ones();
    Point noGhost = Point::Zeros();
    phid.define(dbl, ghostPt);
    phif.define(dbl, ghostPt);
    lapd.define(dbl, noGhost);
    lapf.define(dbl, noGhost);

  }
  {
    PR_TIME("SINGLE_precision_laplacian");
    applyLaplacians<float >(phif, lapf, dbl, numapplies);
  }

  {
    PR_TIME("DOUBLE_precision_laplacian");
    applyLaplacians<float >(phid, lapd, dbl, numapplies);
  }


  PR_TIMER_REPORT();

}  
