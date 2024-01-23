//=========================================
// Dummy implementations
//=========================================

namespace Proto
{
  template<CENTERING Center, typename T, unsigned int C, unsigned int D, unsigned int E>
  EBData<Center,T,C,D,E>::EBData(const Box& a_box, const Irreg<Center>& a_irreg)
  :BoxData<T,C,D,E>{a_box},m_geom{a_irreg} { }
  
template<CENTERING Center, typename T, unsigned int C, unsigned int D, unsigned int E>
EBData<Center,T,C,D,E>::EBData(const Box& a_box, shared_ptr<GeometryService>& a_embeddedBoundary)
  :BoxData<T,C,D,E>{a_box},m_geom{a_embeddedBoundary->irreg<Center>(a_box)} { }


  std::array<double,DIM>  AllRegular::moment(const EBIndex& a_index, int order) const
  {
    static std::array<double, DIM> junk;
    return junk;
  }
  
  std::array<double,DIM>  AllRegular::normal(const EBIndex& a_index) const
  {
    static std::array<double, DIM> junk;
    return junk;
  }

   
  template<typename T, CENTERING DCENTER, CENTERING RCENTER>
  EBDictionary<T,DCENTER,RCENTER>::EBDictionary(std::shared_ptr<GeometryService> a_geom):m_geom(a_geom) {;}

  template<typename T, CENTERING DCENTER, CENTERING RCENTER>
  EBDictionary<T,DCENTER,RCENTER>::EBDictionary(const EBDictionary& a_prior, std::shared_ptr<GeometryService> a_geom):m_geom(a_geom) {;}

  template<typename T, CENTERING DCENTER, CENTERING RCENTER>
  void EBDictionary<T,DCENTER,RCENTER>::registration(const std::string& name,
                      const Stencil<T>&  a_regStencil,
                      const Discretizer<T, DCENTER,RCENTER>& a_ebapprox,
                      const Box& a_span)
  {
    m_map.emplace(name, std::make_tuple(a_ebapprox, a_regStencil, a_span));  
  }

  template<typename T, CENTERING DCENTER, CENTERING RCENTER>
  EBStencil<T,DCENTER,RCENTER> EBDictionary<T, DCENTER, RCENTER>::stencil(const StencilKey& key) const
  {
    auto elem = m_map.at(key.key.first);
    const Box& b = std::get<2>(elem);
    return EBStencil<T, DCENTER, RCENTER>(key.key.second,std::get<1>(elem), b, m_geom,std::get<0>(elem)); }


  // If Stencil operations are always contractive then the input *might* have enough information
  // to construct it's output data holder.  Or, the EBStencil knows the data layout of it's RCENTER
  // object
  template<typename T, CENTERING DCENTER, CENTERING RCENTER>
  EBData<RCENTER,T,1> EBStencil<T, DCENTER, RCENTER>::operator()(const EBData<DCENTER,T,1>& a_src,
                                                                 bool a_initToZero,
                                                                 const T a_scale) const
  {
    return EBData<RCENTER, T, 1>(m_range, m_irreg);
  }

  template<CENTERING Center>
  Irreg<Center>::Irreg() {   }

  template<typename T, CENTERING DCENTER, CENTERING RCENTER>
  EBStencil<T, DCENTER, RCENTER>::EBStencil(const Box& a_rangeBox,
                                            const Stencil<T>& a_regStencil,
                                            const Box& a_span,
                                            const std::shared_ptr<GeometryService>& a_moments,
                                            const Discretizer<T, DCENTER, RCENTER>& a_visitor)
    :m_range(a_rangeBox), m_irreg(a_moments->irreg<RCENTER>(a_rangeBox))
  {
  }
    
}
  
