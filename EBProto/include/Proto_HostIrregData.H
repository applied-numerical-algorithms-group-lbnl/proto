
#ifndef _HostIrregData_H_
#define _HostIrregData_H_

using std::vector;
using std::shared_ptr;

namespace Proto
{
///
/**
*/
  template <class data_t, class index_t, unsigned int ncomp>
  class HostIrregData 
  {
  public:
    ///
    HostIrregData()
    {
    }

    ///
    HostIrregData( const vector<index_t> & a_indicies)
    {
      define(a_indicies);
    }

    ///
    inline void define(const shared_ptr<vector<index_t> >& a_indicies)
    {
      m_indicies = a_indicies;
      std::size_t count = this->size();

      m_data=::std::shared_ptr<T>(new T [count], [](T* p) { delete[] p;});
      m_rawPtr = m_data.get();
    }

    ///
    inline std::size_t size() const {return m_indicies->size()*C;};

    ///
    shared_ptr< vector<index_t> > getIndicies() const
    {
      return m_indicies;
    }


    ///
    T* data()
    {
      return m_rawPtr;
    }

    ///
    const T* data() const
    {
      return m_rawPtr;
    }

    ...
    /**
     */
    virtual data_t* getIndex(const index_t& a_vof,const int& a_comp) const
    {
      PR_assert(a_comp < ncomp)
      data_t* dataPtr =  null_ptr;
      bool found      = false;
      data_t ptrloc   = m_rawPtr;
      for (unsigned int i=0; i<m_indicies->size(); ++i)
      {
        if (a_vof == (*m_indicies)[i]) 
        {
          found = true;
          dataPtr = ptrloc;
          break;
        }
        ptrloc++;
      }

      PR_ASSERT(found, "index not found for hostirregdata");

      if(found)
      {
        dataPtr = ptrloc + a_comp*m_vofs->size();
      }
      return dataPtr;
    }


    ///
    T& operator() (const index_t& a_vof, int a_comp)
    {  return *(getIndex(a_vof, a_comp));}

    ///
    const T& operator() (const index_t& a_vof, int a_comp) const
    {
      return *(getIndex(a_vof, a_comp));
    }

    ///
    inline void setVal(const data_t& a_val)
    {
      for(int ival = 0; ival < size(); ival++)
      {
        m_rawPtr[ival] = a_val;
      }
    }

  protected:

    shared_ptr<vector<index_t> >  m_indicies;
    shared_ptr<data_t>            m_data;
    data_t*                       m_rawPtr;
  };
}

#endif
