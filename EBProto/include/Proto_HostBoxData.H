#ifndef _PROTO_HostBoxData_H_
#define _PROTO_HostBoxData_H_

#include "Proto_Box.H"

namespace Proto
{
/// Multidimensional Rectangular Array that lives on the host
/**
   this is mainly for graph-like data thave lives on the host
*/
  template <class data_t=double, unsigned int ncomp=1>
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
    void copy(const HostBoxData<data_t,Csrc>  & a_dsrc,
              const Box                  & a_srcBox,
              unsigned int                 a_srcComp,
              const Box                  & a_destBox,
              unsigned int                 a_destComp,
              unsigned int                 a_numcomp);
    
    ///this DOES include sizeof(data_t)  
    virtual std::size_t charsize(const Box&   a_bx,
                                 unsigned int a_startcomp,
                                 unsigned int a_numcomps) const
    {
      PR_assert((a_startcomp + a_numcomps) < ncomp);
      Box binter = a_bx & m_box;
      size_t  numpts = binter.size();
      
      size_t regsize = numpts*a_numcomps*sizeof(data_t);

      return regsize;
      
    }
    virtual void  linearIn(void*                  a_buf,
                           const Box       &      a_bx,
                           unsigned int           a_startcomp,
                           unsigned int           a_numcomps)
    {
      PR_assert((a_startcomp + a_numcomps) < ncomp);

      char* charbuf = (char*)a_buf;
      Box binter = a_bx & m_box;
      for(auto bit = binter.begin(); bit != binter.end(); ++bit)
      {
        Point pt = *bit;
        for(int icomp = 0; icomp < a_numcomps; icomp++)
        {
          data_t* dataptr = (data_t*)  charbuf;
          (*this)(pt, icomp + a_startcomp) = *dataptr;
          charbuf += sizeof(data_t);
        }
      }


    }
    virtual void linearOut(void*               a_buf,
                           const Box       &   a_bx,
                           unsigned int        a_startcomp,
                           unsigned int        a_numcomps) const
    {
      PR_assert((a_startcomp + a_numcomps) < ncomp);

      char* charbuf = (char*)a_buf;
      Box binter = a_bx & m_box;
      for(auto bit = binter.begin(); bit != binter.end(); ++bit)
      {
        Point pt = *bit;
        for(int icomp = 0; icomp < a_numcomps; icomp++)
        {
          data_t* dataptr = (data_t*)  charbuf;
          *dataptr = (*this)(pt, icomp + a_startcomp);
          charbuf += sizeof(data_t);
        }
      }
    }

    ///
    inline const data_t& operator()(const Point& a_pt, 
                               unsigned int  a_c) const;

    ///
    inline data_t& operator()(const Point& a_pt, 
                         unsigned int  a_c);
  
    ///
    inline Box box() const {return m_box;};
    
    ///
    inline std::size_t size() const {return m_box.size()*ncomp;};

    ///
    inline bool defined() const {return bool(m_data);};
   

    ///
    void setVal(const data_t& a_val);
  

    ///
    data_t* data()
    {
      return m_rawPtr;
    }

    ///
    const data_t* data() const
    {
      return m_rawPtr;
    }

  private:
    /// Copy constructor/assignment forbidden for all the usual reasons.
    HostBoxData(const HostBoxData<data_t,ncomp>& a_src);
    HostBoxData& operator=(const HostBoxData<data_t,ncomp>& a_src);
    
    Box                       m_box;    ///< Box defining the domain of *this
    ::std::shared_ptr<data_t> m_data;   ///< Data array
    data_t*                   m_rawPtr; ///< Raw pointer to the data
  };

} //end Proto namespace

#include "Proto_HostBoxDataImplem.H"
#endif //end include guard
