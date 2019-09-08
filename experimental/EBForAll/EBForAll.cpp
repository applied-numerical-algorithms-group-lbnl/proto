
#include "Proto.H"
#include "EBProto.H"

using namespace Proto;



template<typename T>
inline T&
getIrregData(T& a_s)
{
  return a_s;
}

template <CENTERING cent, typename  data_t, unsigned int ncomp>
inline IrregData<cent, data_t, ncomp>&
getIrregData(EBBoxData<cent, data_t, ncomp>& a_s)
{
  return a_s.getIrregData();
}


template< typename data_t, unsigned int ncomp>
struct
uglyStruct
{
  data_t*      m_startPtr;
  unsigned int m_varsize;
  unsigned int m_index;
};

template <CENTERING cent, typename T>
inline T
getUglyStruct(const vector<EBIndex<cent> >& a_indices,
              T& a_T)
{
  return a_T;
}

template <CENTERING cent, typename  data_t, unsigned int ncomp>
inline vector< uglyStruct<data_t, ncomp> >
getUglyStruct(const vector<EBIndex<cent> >& a_indices,
              IrregData<cent, data_t, ncomp>& a_s )
{
  vector< uglyStruct<data_t, ncomp> > retval;
  for(int ivec = 0; ivec < a_indices.size(); ivec++)
  {
    uglyStruct<data_t, ncomp>  vecval;
    vecval.m_startPtr = a_s.data();
    vecval.m_varsize  = a_s.vecsize();
    vecval.m_index    = a_s.index(a_indices[ivec], 0);
    retval.push_back(vecval);
  }
  return retval;
}


#ifdef PROTO_CUDA
template<typename Func, typename... Srcs>
inline void
cudaEBforAll(const Func & a_F,  Box a_box, Srcs&... a_srcs)
{
}
#else
///
template<typename T>
inline T
getVar(unsigned int ivec,  T a_s)
{
  return a_s;
}
///
template<typename data_t, unsigned int ncomp>
inline Var<data_t, ncomp>
getVar(unsigned int a_ivec,
       vector< uglyStruct<data_t, ncomp> > a_dst)
{
  Var<data_t, ncomp> retval;
  const uglyStruct<data_t, ncomp> ugly = a_dst[a_ivec];
  for(int icomp = 0; icomp < ncomp; icomp++)
  {
    retval.m_ptrs[icomp] = ugly.m_startPtr + ugly.m_index + (ugly.m_varsize*icomp);
  }
  return retval;
}
///going into this srcs are vector<uglystruct> and other stuff
template<typename data_t,unsigned int ncomp,  typename Func, typename... Srcs>
void
vectorFunc(const Func& a_F, vector< uglyStruct<data_t, ncomp> > a_dst, Srcs... a_srcs)
{
  for(unsigned int ivec = 0; ivec < a_dst.size(); ivec++)
  {
    a_F(getVar(ivec, a_dst), (getVar(ivec, a_srcs))...);
  }
       
}

///going into this srcs are IrregDatas and other stuff
template<CENTERING cent, typename  data_t, unsigned int ncomp, typename Func, typename... Srcs>
inline void
hostEBForAllIrreg(const Func& a_F, const Box& a_box,
                  IrregData<cent, data_t, ncomp>& a_dst,
                  Srcs&...  a_srcs)
{
  //indicies into irreg vector that correspond to input box
  vector<EBIndex<cent> > dstvofs = a_dst.getIndices(a_box);
  vectorFunc(a_F, getUglyStruct(dstvofs, a_dst), (getUglyStruct(dstvofs, a_srcs))...);
}

template<typename T>
inline T&
getBoxData(T& a_s)
{
  return a_s;
}

template <CENTERING cent, typename  data_t, unsigned int ncomp>
inline BoxData<data_t, ncomp>&
getBoxData(EBBoxData<cent, data_t, ncomp>& a_s)
{
  return a_s.getRegData();
}

///going into this srcs are EBBoxDatas and other stuff
template<typename Func, typename... Srcs>
inline void
hostEBforAll(const Func & a_F,  Box a_box, Srcs&... a_srcs)
{
  //call regular forall
  forallInPlaceBase(a_F, a_box, (getBoxData(a_srcs))...);
  
  //do the same thing for the irregular data
  hostEBForAllIrreg(a_F, a_box, getIrregData(a_srcs)...);
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

typedef Var<double,DIM> V;


PROTO_KERNEL_START 
unsigned int setUF(V a_U, double  a_val)
{
  for(int idir = 0; idir < DIM; idir++)
  {
    a_U(idir) = a_val;
  }
}
PROTO_KERNEL_END(setUF, setU)


PROTO_KERNEL_START 
unsigned int setVF(V a_V, double  a_val)
{
  for(int idir = 0; idir < DIM; idir++)
  {
    a_V(idir) = a_val;
  }
}
PROTO_KERNEL_END(setVF, setV)


PROTO_KERNEL_START 
unsigned int setWtoUplusVF(V a_W,
                           V a_U,
                           V a_V,
                           double  a_val)
{
  for(int idir = 0; idir < DIM; idir++)
  {
    a_W(idir) = a_U(idir) + a_V(idir);
    if(a_W(idir) != a_val)
    {
      printf("values do not match");
    }
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
  shared_ptr<BaseIF>              impfunc(new SimpleEllipsoidIF(ABC, X0, R, false));
  shared_ptr<GeometryService<2> > geoserv(new GeometryService<2>(impfunc, origin, dx, domain, grids, geomGhost, 0));
  shared_ptr<LevelData<EBGraph> > graphs = geoserv->getGraphs(domain);

  for(int ibox = 0; ibox < grids.size(); ibox++)
  {
    Box grid = grids[ibox];
    EBBoxData<CELL, double, DIM> U(grid,(*graphs)[ibox]);
    EBBoxData<CELL, double, DIM> V(grid,(*graphs)[ibox]);
    EBBoxData<CELL, double, DIM> W(grid,(*graphs)[ibox]);
    unsigned long long int numFlopsPt = 0;
    double uval = 1;
    double vval = 2;
    double wval = 3;
    EBforAllInPlace(numFlopsPt, "setU", setU, grid, U, uval);
    EBforAllInPlace(numFlopsPt, "setV", setV, grid, V, vval);
    EBforAllInPlace(numFlopsPt, "setWtoUPlusV", setWtoUplusV, grid, W, U, V, wval);


  }

}








