#ifndef _PROTO_HostBoxData_H_
#define _PROTO_HostBoxData_H_

#include "Proto_Box.H"

namespace Proto
{
/// Multidimensional Rectangular Array that lives on the host
/**
   this is mainly for graph-like data thave lives on the host
*/
  template <class T=double, unsigned int C=1>
  class HostBoxData 
  {
  public:

    ///
    HostBoxData();

    ///
    HostBoxData(const Box& a_box);
    
    /// 
    void define(const Box& a_box);
   
    /// Destructor.
    ~HostBoxData(); 

    ///
    template<unsigned int Csrc>
    void copy(const HostBoxData<T,Csrc>  & a_dsrc,
              const Box                  & a_srcBox,
              unsigned int                 a_srcComp,
              const Box                  & a_destBox,
              unsigned int                 a_destComp,
              unsigned int                 a_numcomp);
    
    ///
    inline const T& operator()(const Point& a_pt, 
                               unsigned int  a_c) const;

    ///
    inline T& operator()(const Point& a_pt, 
                         unsigned int  a_c);
  
    ///
    inline Box box() const {return m_box;};
    
    ///
    inline std::size_t size() const {return m_box.size()*C;};

    ///
    inline bool defined() const {return bool(m_data);};
   

    ///
    void setVal(const T& a_val);
  

  private:
    
    Box                     m_box;    ///< Box defining the domain of *this
    ::std::shared_ptr<T>    m_data;   ///< Data array
    T*                      m_rawPtr; ///< Raw pointer to the data
  };

} //end Proto namespace

#include "Proto_HostBoxDataImplem.H"
#endif //end include guard
