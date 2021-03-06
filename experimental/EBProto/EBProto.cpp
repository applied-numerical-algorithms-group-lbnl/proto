
#include "Proto.H"
#include "Proto_EBProto.H"

/*
namespace Proto
{
  struct EBIndex
  {
    Point index;
    unsigned int  i;
  };

  enum CENTERING {CELL,XFACE,YFACE,ZFACE};

  class GeometryService;  //forward declaration
  
  template<CENTERING Center>
  struct Irreg
  {
    Irreg(Irreg&&)=default;
    Irreg(const Irreg& a_irreg)=default;
    
  private:
    friend class GeometryService;
    Irreg();
  };
  
  class GeometryService
   {
   public:
     virtual ~GeometryService()=default;
     virtual std::array<double,DIM> moment(const EBIndex& a_index, int order) const=0;
     virtual std::array<double,DIM> normal(const EBIndex& a_index) const=0;
     template<CENTERING Center>
     Irreg<Center> irreg(const Box& a_box) {return Irreg<Center>();}

     // still need a solution to the templated virtual function problem (bvs)
     //virtual Irreg<Center> getIrreg(const Box& a_box) {return Irreg<Center>(); }
   };

  class AllRegular : public GeometryService
  {
  public:
    virtual ~AllRegular()=default;
    virtual std::array<double,DIM> moment(const EBIndex& a_index, int order) const;
    virtual std::array<double,DIM> normal(const EBIndex& a_index) const;
  };
  


 
  
  template<CENTERING Center=CELL, typename T=double, unsigned int C=1, unsigned char D=0, unsigned char E=0>
  class EBData : private BoxData<T,C,D,E>
  {
  public:
    // using  center_type = Center;
    using  value_type  = T;
    EBData(const Box& a_box, shared_ptr<GeometryService>& a_embeddedBoundary);
    EBData(const Box& a_box, const Irreg<Center>& a_irreg);
    EBData(EBData&&)=default;

    // aliasing and member functions replicated from BoxData
  private:
    Irreg<Center> m_geom;
  };

  template<typename T, unsigned char DType, unsigned char RType>
  class Discretizer
  {
  public:
    virtual ~Discretizer()=default;
    Discretizer()=default;
    // returns single-valued and multi-valued operators
  //   virtual std::pair<Stencil<T>, std::vector<std::pair<EBIndex,T>>> irregular(const EBIndex& a_point,
  //                                                                             const GeometryService& a_moments,
  //                                                                             const std::set<EBIndex>& a_reachable,
  //                                                                             const Box& a_span);
  //
  };
    
  template<typename T, CENTERING DCENTER, CENTERING RCENTER>
  class EBStencil
  {
  public:
    /// strong constructor
    EBStencil(const Box& a_rangeBox, const Stencil<T>& a_regStencil, const Box& a_span,
              const std::shared_ptr<GeometryService>& a_moments,
              const Discretizer<T, DCENTER, RCENTER>& a_visitor);
    
    EBStencil(EBStencil&&) = default;
    
    /// remaining construction takes C++ defaults
    
    void apply(const EBData<DCENTER,T, 1>&   a_src,
               EBData<RCENTER,T, 1>&         a_dest,
               const Box&     a_bx,
               bool           a_initToZero = false,
               const T        a_scale = 1) const { ; }

    // Box-inference / move semantics means that EBStencil can create other CENTERINGS
    EBData<RCENTER,T,1> operator()(const EBData<DCENTER,T,1>& a_src,
                                   bool a_initToZero = false,
                                   const T a_scale = 0) const ; // full box inference version
    
    // stencil algebra elided for clarity.  Recapitulates Stencil
  private:
    Box m_range;
    Irreg<RCENTER> m_irreg;
  };


  // example class of Discretizer
  template<typename T, CENTERING DCENTER, CENTERING RCENTER>
  class QuadraticPoly: public Discretizer<T, DCENTER, RCENTER>
  {
  public:
    QuadraticPoly()=default;
    // virtual std::pair<Stencil<T>, std::vector<std::pair<EBIndex,T>>> irregular(const EBIndex& a_point,
    //                                                                           const GeometryService& a_moments,
    //                                                                           const std::set<EBIndex>& a_reachable,
    //                                                                           const Box& a_span);
  };

  struct StencilKey
  {
    std::pair<std::string,Box> key;
    StencilKey(std::string skey,Box bkey)
    {
      key=std::make_pair(skey,bkey);
    }
  };



  template<typename T, CENTERING DCENTER, CENTERING RCENTER>
  class EBDictionary
  {
  public:
    EBDictionary(std::shared_ptr<GeometryService> a_geom):m_geom(a_geom) {;}
    EBDictionary(const EBDictionary& a_prior, std::shared_ptr<GeometryService> a_newGeom)
{;
  //migrate
}

    void registration(const std::string& name,
                      const Stencil<T>&  a_regStencil,
                      const Discretizer<T, DCENTER,RCENTER>& a_ebapprox,
                      const Box& a_span)
    {
      m_map.emplace(name, std::make_tuple(a_ebapprox, a_regStencil, a_span));
      //auto m = m_map[name];
      //std::get<0>(m)=a_ebapprox;
      //m_map[name]=std::make_tuple(a_ebapprox, a_regStencil, a_span);
    }
    

    EBStencil<T,DCENTER,RCENTER> stencil(const StencilKey& key) const;
  

  private:
    std::shared_ptr<GeometryService> m_geom;
    std::map<std::string,std::tuple<const Discretizer<T,DCENTER,RCENTER>&, const Stencil<T>&, const Box >> m_map;
    
  };
  
}

*/
using namespace Proto;

int main(int argc, char* argv[])
{

  std::shared_ptr<GeometryService> geom = std::shared_ptr<GeometryService>(new AllRegular());
  EBDictionary<double, CELL, XFACE> eb_xf(geom);
  
  Stencil<double> average(Shift::Zeros(),1.0);
  average += 1.0*Shift::Basis(0,-1);
  average*=0.5;

  eb_xf.registration("Average", average,                          
                     QuadraticPoly<double, CELL, XFACE>(),
                     average.span().grow(1));
  Point lo(-12,-15,-10,-14);
  Point hi(34,12,6,12);
  EBIndex eb1{lo,0};
  EBIndex eb2{Point(4,3,6,7,8),1};
  Box d{lo, hi};
  Box r{lo+Point::Basis(0), hi};

  auto stencil = eb_xf.stencil({"Average",r});
  EBData<CELL>  cellData(d,geom);
  auto xface = stencil(cellData); //yeah bitches, inference for EBData objects.

}








