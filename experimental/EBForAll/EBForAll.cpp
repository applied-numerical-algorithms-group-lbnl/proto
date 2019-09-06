
#include "Proto.H"
#include "EBProto.H"

using namespace Proto;
template<typename data_t, unsigned int C>
class EBVar
{
public:
  data_t vars[C];
  data_t& operator() (int a_index)
  {
    return vars[a_index];
  }
  const data_t& operator() (int a_index) const
  {
    return vars[a_index];
  }
};


template<typename Func, typename... Srcs>
inline void hostEBforAll(const Func & a_F,  Box a_box, Srcs&... a_srcs)
{
}


#ifdef PROTO_CUDA
template<typename Func, typename... Srcs>
inline void hostEBforAll(const Func & a_F,  Box a_box, Srcs&... a_srcs)
{
}
#endif

template<typename Func, typename... Srcs>
inline void EBforAllInPlace(unsigned long long int a_num_flops_point,
                            const char*            a_timername,
                            const Func & a_F,  Box a_box, Srcs&... a_srcs)
{
  PR_TIME(a_timername);


  unsigned long long int boxfloops = a_num_flops_point*a_box.size();

#ifdef PROTO_CUDA
  cudaEBforAll(a_F, a_box, a_srcs...);
#else
  hostEBforAll(a_F, a_box, a_srcs...);
#endif
  PR_FLOPS(boxfloops);
}

typedef EBVar<double,DIM> V;


PROTO_KERNEL_START 
unsigned int setUF(V& a_U)
{
  for(int idir = 0; idir < DIM; idir++)
  {
    a_U(idir) = idir;
  }
}
PROTO_KERNEL_END(setUF, setU)


PROTO_KERNEL_START 
unsigned int setVF(V& a_V)
{
  for(int idir = 0; idir < DIM; idir++)
  {
    a_V(idir) = 2*idir;
  }
}
PROTO_KERNEL_END(setVF, setV)


PROTO_KERNEL_START 
unsigned int setWtoUplusVF(V& a_W,
                           const V& a_U,
                           const V& a_V)
{
  for(int idir = 0; idir < DIM; idir++)
  {
    a_W(idir) = a_U(idir) + a_V(idir);
  }
}
PROTO_KERNEL_END(setWtoUplusVF, setWtoUplusV)

int main(int argc, char* argv[])
{

  int nx = 16;
  int maxGrid = nx; //this test is for a single grid
  Box domain(Point::Zeros(), Point::Ones(nx-1));
  double dx = 1.0/domain.size(0);
  std::array<bool, DIM> periodic;
  for(int idir = 0; idir < DIM; idir++) periodic[idir]=true;
  DisjointBoxLayout grids(domain, maxGrid, periodic);

  Point dataGhost = Point::Zeros();
  Point geomGhost = Point::Zeros();

  RealVect ABC = RealVect::Unit();
  RealVect X0 = RealVect::Unit();
  X0 *= 0.5;
  RealVect origin= RealVect::Zero();
  double R = 0.25;
  shared_ptr<BaseIF>                       impfunc(new SimpleEllipsoidIF(ABC, X0, R, false));
  shared_ptr<GeometryService<2> >  geoserv(new GeometryService<2>(impfunc, origin, dx, domain, grids, geomGhost, 0));
  shared_ptr<LevelData<EBGraph> > graphs = geoserv->getGraphs(domain);

  for(int ibox = 0; ibox < grids.size(); ibox++)
  {
    Box grid = grids[ibox];
    EBBoxData<CELL, double, DIM> U(grid);
    EBBoxData<CELL, double, DIM> V(grid);
    EBBoxData<CELL, double, DIM> W(grid);
    unsigned long long int numFlopsPt = 0;
    EBforAllInPlace(numFlopsPt, "setU", setU, grid, U);
    EBforAllInPlace(numFlopsPt, "setV", setU, grid, U);
    EBforAllInPlace(numFlopsPt, "setWtoUPlusV", setU, grid, W, U, V);


  }

}








