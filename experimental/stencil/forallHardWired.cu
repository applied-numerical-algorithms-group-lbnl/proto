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
  a_nmult= 20;
  cout << "nx          = " << a_nx << endl;
  cout << "num_applies = " << a_numapplies << endl; 
  cout << "maxgrid     = " << a_maxgrid << endl;
  cout << "num_mult is hardwired to 20" << endl;
}


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
__global__
void hardwiredRiemann(unsigned int a_Nz, unsigned int a_zinc, unsigned int a_varinc,
                      double* a_out, double* a_low, double* a_hig, unsigned int a_idir, double gamma, int a_nmult)
{
  unsigned int idx = threadIdx.x + blockIdx.x*blockDim.x;
  double* out = a_out + idx;
  double* hig = a_hig + idx;
  double* low = a_low + idx;
  unsigned int roff = 0;
  unsigned int uoff = a_varinc*(a_idir + 1);
  unsigned int poff = a_varinc*(NUMCOMPS-1);
  for(unsigned int zloc = 0; zloc < a_Nz; zloc++)
  {

    double& rhoo = *(out + roff);
    double& rhol = *(low + roff);
    double& rhor = *(hig + roff);
    double& uo   = *(out + uoff);
    double& ul =   *(low + uoff);
    double& ur =   *(hig + uoff);
    double& po   = *(out + poff);
    double& pl =   *(low + poff);
    double& pr =   *(hig + poff);


    //2
    double rhobar = (rhol + rhor)*.5;
    //2
    double pbar = (pl + pr)*.5;
    //2
    double ubar = (ul + ur)*.5;
    //took this one out for a bunch of multiplies so
    //I can have flops I can count
//  double cbar = sqrt(gamma*pbar/rhobar);
    //2
    double cbar = gamma*pbar/rhobar;
    //NMULT
    cbar *= pbar;
    cbar *= pbar;
    cbar *= pbar;
    cbar *= pbar;
    cbar *= pbar;
    cbar *= pbar;
    cbar *= pbar;
    cbar *= pbar;
    cbar *= pbar;
    cbar *= pbar;
    cbar *= pbar;
    cbar *= pbar;
    cbar *= pbar;
    cbar *= pbar;
    cbar *= pbar;
    cbar *= pbar;
    cbar *= pbar;
    cbar *= pbar;
    cbar *= pbar;
    cbar *= pbar;

    //7
    double pstar = (pl + pr)*.5 + rhobar*cbar*(ul - ur)*.5;
    //7
    double ustar = (ubar + ur)*.5 + (pl - pr)/(2*rhobar*cbar);

    rhoo  = rhol;
    uo    =   ul;
    po    =   pl;
    //4
    rhoo += (pstar - po)/(cbar*cbar);
    uo    = ustar;
    po = pstar;


    out += a_zinc;
    hig += a_zinc;
    low += a_zinc;
  }
}

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
    PR_TIME("hardwired_riemann");
    cout << "doing riemann problems " << endl;
    for(unsigned int ibox = 0; ibox < localBoxes.size(); ibox++)
    {
      for(unsigned int iapp = 0; iapp < a_numapplies; iapp++)
      {
        int istream = iapp%a_numstream;
        Box appBox       = a_dbl[localBoxes[ibox]];
        
        int stride = appBox.size(0);
        int blocks = appBox.size(1);
         
#if DIM==3
        unsigned int Nz        = appBox.size(2);
        unsigned int zinc      = appBox.flatten(2).size();
        unsigned int varinc    = appBox.size(); //has to be non-zero or we have an infinite loop
#else

        unsigned int Nz        = 1;
        unsigned int zinc      = 1; //has to be non-zero or we have an infinite loop
        unsigned int varinc    = appBox.size(); //has to be non-zero or we have an infinite loop
#endif
        unsigned long long int count = (28 + a_nmult)*appBox.size();
        size_t smem = 0;
        hardwiredRiemann<<<blocks, stride, smem, streams[istream]>>>(Nz, zinc, varinc, a_out[ibox].data(), a_low[ibox].data(), a_hig[ibox].data(), idir, gamma, a_nmult);

        PR_FLOPS(count);
      }
    }
    sync();
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
    PR_TIME("1_STREAM");
    int nstream = 1;
    cout << "running test with " << nstream << " stream(s)" << endl;
    doSomeForAlls(out, hig, low, dbl, niter, nstream, nmult);
  }
  {
    PR_TIME("2_STREAMS");
    int nstream = 2;
    cout << "running test with " << nstream << " stream(s)" << endl;
    doSomeForAlls(out, hig, low, dbl, niter, nstream, nmult);
  }

  {
    PR_TIME("4_STREAMS");
    int nstream = 4;
    cout << "running test with " << nstream << " stream(s)" << endl;
    doSomeForAlls(out, hig, low, dbl, niter, nstream, nmult);
  }

  {
    PR_TIME("8_STREAMS");
    int nstream = 8;
    cout << "running test with " << nstream << " stream(s)" << endl;
    doSomeForAlls(out, hig, low, dbl, niter, nstream, nmult);
  }


  {
    PR_TIME("16_STREAMS");
    int nstream = 16;
    cout << "running test with " << nstream << " stream(s)" << endl;
    doSomeForAlls(out, hig, low, dbl, niter, nstream, nmult);
  }

  {
    PR_TIME("32_STREAMS");
    int nstream = 32;
    cout << "running test with " << nstream << " stream(s)" << endl;
    doSomeForAlls(out, hig, low, dbl, niter, nstream, nmult);
  } 


  PR_TIMER_REPORT();

}  
