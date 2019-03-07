
#include "Proto.H"
#include <array>
#include <utility>
#include <set>

namespace Proto
{
  struct PIndex
  {
    Point index;
    unsigned int   vol;
  };

  class GeometryService
   {
   public:
     virtual std::array<double,DIM> moment(const PIndex& a_index, int order) const=0;
     virtual std::array<double,DIM> normal(const PIndex& a_index) const=0;
   };

  class AllRegular : public GeometryService
  {
  public:
    virtual std::array<double,DIM> moment(const PIndex& a_index, int order) const;
    virtual std::array<double,DIM> normal(const PIndex& a_index) const;
  };
  
  enum CENTERING {CELL,XFACE,YFACE,ZFACE};
    
  template<CENTERING Center=CELL, typename T=double, unsigned int C=1, unsigned char D=0, unsigned char E=0>
  class EBData : private BoxData<T,C,D,E>
  {
  public:
    // using  center_type = Center;
    using  value_type  = T;
    EBData(const Box& a_box, const GeometryService& a_embeddedBoundary);

    // aliasing and member functions replicated from BoxData
  };

  template<typename T, unsigned char DType, unsigned char RType>
  class Discretizer
  {
  public:
    // returns single-valued and multi-valued operators
    virtual std::pair<Stencil<T>, std::vector<std::pair<PIndex,T>>> irregular(const PIndex& a_point,
                                                                              const GeometryService& a_moments,
                                                                              const std::set<PIndex>& a_reachable,
                                                                              const Box& a_span);
  };
    
  template<typename T, CENTERING DCENTER, CENTERING RCENTER>
  class EBStencil
  {
  public:
    /// strong constructor
    EBStencil(const Box& a_rangeBox, const Stencil<T>& a_regStencil, const Box& a_span,
              const GeometryService& a_moments,
              const Discretizer<T, DCENTER, RCENTER>& a_visitor);
    
    /// remaining construction takes C++ defaults
    
    void apply(const EBData<DCENTER,T, 1>&   a_src,
               EBData<RCENTER,T, 1>&         a_dest,
               const Box&     a_bx,
               bool           a_initToZero = false,
               const T        a_scale = 1) const;

    // stencil algebra elided for clarity.  Recapitulates Stencil
  };


  // example class of Discretizer
  template<typename T, CENTERING DCENTER, CENTERING RCENTER>
  class QuadraticPoly: public Discretizer<T, DCENTER, RCENTER>
  {
    virtual std::pair<Stencil<T>, std::vector<std::pair<PIndex,T>>> irregular(const PIndex& a_point,
                                                                              const GeometryService& a_moments,
                                                                              const std::set<PIndex>& a_reachable,
                                                                              const Box& a_span);
  };
}

using namespace Proto;

int main(int argc, char* argv[])
{

  AllRegular geom;
  Stencil<double> average(Shift::Zeros(),1.0);
  average += 1.0*Shift::Basis(0,-1);
  
  average*=0.5;
  Point lo(-12,-15,-10,-14);
  Point hi(34,12,6,12);
  Box d{lo, hi};
  Box r{lo+Point::Basis(0), hi}; //can't deduce this one from what I can tell

  EBStencil<double, CELL, XFACE> stencil(r, average, average.span().grow(1), geom,
                                         QuadraticPoly<double, CELL, XFACE>());
  
  EBData<CELL>  cellData(d,geom);
  EBData<XFACE> xface(r,geom);
  stencil.apply(cellData, xface, d);

}
  


