#ifndef _Proto_REGISTER_H_
#define _Proto_REGISTER_H_
#include "Proto.H"
#include "Proto_MemType.H"
namespace Proto
{
  template<typename T, unsigned int C, MemType MEMTYPE>
  struct Register
  {
    shared_ptr<BoxData<T,C,MEMTYPE> > m_data;
    unsigned int m_dir;
    unsigned int m_lohi;
    Register(
             shared_ptr<BoxData<T,C,MEMTYPE> > a_data,
             unsigned int a_dir,
             unsigned int a_lohi)
    {
      m_data = a_data;
      m_dir = a_dir;
      m_lohi = a_lohi;
    };
  };
}
#endif
