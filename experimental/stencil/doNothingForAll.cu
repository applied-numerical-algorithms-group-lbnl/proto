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
using std::cout;
using std::endl;
using std::vector;
using std::shared_ptr;
using namespace Proto;
constexpr unsigned int NUMCOMPS=DIM+2;
typedef Var<double,NUMCOMPS> State;

//arbitrary int for number of multiplies to replace square root in riemann

/**/
void
parseCommandLine(unsigned int& a_nmult, unsigned int & a_nx, unsigned int& a_maxgrid, unsigned int & a_numapplies, int argc, char* argv[])
{
  //defaults
  a_nx = 256;
  a_numapplies = 100;
  a_nmult = 100;
  a_maxgrid = 64;
  cout << "kernel timings of various foralls ([] shows defaults)" << endl;
  cout << "usage:  " << argv[0] << " -n nx[256] -m max_grid[128] -a num_iterations[100] -s max_streams[32] -t num_multiplies[100]" << endl;
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
    else if(strcmp(argv[iarg], "-t") == 0)
    {
      a_nmult = atoi(argv[iarg+1]);
    }
  }
  cout << "nx          = " << a_nx << endl;
  cout << "num_applies = " << a_numapplies << endl; 
  cout << "maxgrid     = " << a_maxgrid << endl;
  cout << "num_mult    = " << a_nmult << endl;
}




PROTO_KERNEL_START
void doNothingF(State& a_out,
                const State& a_low,
                const State& a_high,
                int   a_dir,
                double a_gamma,
                unsigned int a_nmult)
{
}
PROTO_KERNEL_END(doNothingF, doNothing)



struct DoNothingStruct 
{ 
  __device__ void op(State& a_out,
                     const State& a_low,
                     const State& a_high,
                     int   a_dir,
                     double a_gamma,
                     unsigned int a_nmult)
  { 
    return doNothing(a_out, a_low, a_high, a_dir, a_gamma, a_nmult);
  }
};


///proxies for Chombo-style SPMD functions
unsigned int CH_numProc()
{
  return 1;
}

unsigned int CH_procID()
{
  return 0;
}

///boxes that cover a domain box
class DisjointBoxLayout
{
private:
  struct localData
  {
    Box                   m_coarsenedDom;
    vector<unsigned int>  m_procs;
    vector<unsigned int>  m_localBoxes;
    unsigned int          m_maxgrid;
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
  unsigned  int procID(unsigned int a_index) const
  {
    PROTO_ASSERT(m_internals,"trying to access undefined dbl procids");
    Point coarpt = m_internals->m_coarsenedDom[a_index];
    return m_internals->m_procs[a_index];
  }

  ///
  Box operator[](unsigned int a_index) const
  {
    PROTO_ASSERT(m_internals,"trying to access undefined dbl boxes");
    Point coarpt = m_internals->m_coarsenedDom[a_index];
    Box coarBox(coarpt, coarpt);
    Box retval = coarBox.refine(m_internals->m_maxgrid);
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
    PROTO_ASSERT(m_internals,"trying to access undefined dbl local boxes");
    return m_internals->m_localBoxes;
  }

};

///data over a disjointboxlayout with ghost cells
template <class T>
class LevelData
{

public:
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

  void setToZero()
  {
    for(unsigned int ibox = 0; ibox < m_data.size(); ibox++)
    {
      m_data[ibox]->setVal(0);
    }
  }
private: 
  //in parallel, this is be the data on this proc
  vector<shared_ptr<T> >        m_data;
  DisjointBoxLayout             m_grids;
  bool                          m_isDefined;
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

void
doSomeForAlls(  LevelData< BoxData<double, NUMCOMPS> > & a_out,
                LevelData< BoxData<double, NUMCOMPS> > & a_low,
                LevelData< BoxData<double, NUMCOMPS> > & a_hig,
                const DisjointBoxLayout & a_dbl,
                const unsigned int      & a_numapplies,
                const unsigned int      & a_numstream,
                const unsigned int      & a_nmult)
{

  //remember this is just for timings
  vector<unsigned int> localBoxes = a_dbl.localBoxes();
  cout << "local boxes size  = " << localBoxes.size() << endl;
  for(unsigned int ibox = 0; ibox < localBoxes.size(); ibox++)
  {
    a_out[ibox].setVal(1.);
    a_hig[ibox].setVal(1.);
    a_low[ibox].setVal(1.);
  }
  vector<cudaStream_t> streams(a_numstream);
  double gamma = 1.4;
  int idir = 0;
  for(unsigned int ibox = 0; ibox < a_numstream; ibox++)
  {
    cudaStreamCreate(&streams[ibox]);
  }


  {
    PR_TIME("No_Z_increment");
    {
      cout << "doing empty foralls " << endl;
      for(unsigned int ibox = 0; ibox < localBoxes.size(); ibox++)
      {
        for(unsigned int iapp = 0; iapp < a_numapplies; iapp++)
        {
          PR_TIME("do_nothing_on_level_multiStream");
          int istream = iapp%a_numstream;
          Box appBox       = a_dbl[localBoxes[ibox]];

          cudaForallStruct(streams[istream], DoNothingStruct()  , appBox, a_out[ibox], a_low[ibox], a_hig[ibox], idir, gamma, a_nmult);
        }
      }
      sync();
    }

  }

  {
    PR_TIME("Z_increment_Version");
    cout << "z incr version" << endl;
    {
      cout << "doing empty foralls " << endl;
      for(unsigned int ibox = 0; ibox < localBoxes.size(); ibox++)
      {
        for(unsigned int iapp = 0; iapp < a_numapplies; iapp++)
        {
          PR_TIME("do_nothing_on_level_multiStream");
          int istream = iapp%a_numstream;
          Box appBox       = a_dbl[localBoxes[ibox]];

          cudaForallZincStruct(streams[istream], DoNothingStruct()  , appBox, a_out[ibox], a_low[ibox], a_hig[ibox], idir, gamma, a_nmult);
        }
      }
      sync();
    }

  }


  {
    PR_TIME("Empty_Indexer_Version");
    cout << "empty indexer version" << endl;
    {
      cout << "doing empty foralls " << endl;
      for(unsigned int ibox = 0; ibox < localBoxes.size(); ibox++)
      {
        for(unsigned int iapp = 0; iapp < a_numapplies; iapp++)
        {
          PR_TIME("do_nothing_on_level_multiStream");
          int istream = iapp%a_numstream;
          Box appBox       = a_dbl[localBoxes[ibox]];

          cudaForallEmptyIndexer(streams[istream], DoNothingStruct()  , appBox, a_out[ibox].box());
        }
      }
      sync();
    }

  }


  for(unsigned int ibox = 0; ibox < a_numstream; ibox++)
  {
    cudaStreamDestroy(streams[ibox]);
  }
}
/**/
int main(int argc, char* argv[])
{
  //have to do this to get a time table
  PR_TIMER_SETFILE("proto.time.table");
  unsigned int nx, niter, maxgrid,  nmult;
  parseCommandLine(nmult, nx, maxgrid, niter,  argc, argv);

  Point lo = Point::Zeros();
  Point hi = Point::Ones(nx - 1);
  Box domain(lo, hi);
  
  DisjointBoxLayout dbl(domain, maxgrid);
  LevelData<BoxData<double, NUMCOMPS> > out, hig, low;

  {
    
    PR_TIME("data definition");

    out.define(dbl, Point::Zeros());
    hig.define(dbl, Point::Zeros());
    low.define(dbl, Point::Zeros());

  }

  {
    PR_TIME("8_STREAMS");
    int nstream = 8;
    cout << "running test with " << nstream << " stream(s)" << endl;
    doSomeForAlls(out, hig, low, dbl, niter, nstream, nmult);
  }



  PR_TIMER_REPORT();

}  
