#pragma once
#include "PR_Real.h"


template<typename Func, unsigned int C>
inline void forall(const int size, Func func, Var<C> v)
{
 // int tripCount=0;
  Var<C> c = v;
  int i=0;
  for(; i<size-VLENGTH; i+=VLENGTH,++c)
  {
    func(c); // tripCount++;
  }

  int fixup=size-VLENGTH*(size/VLENGTH);
  if(fixup>0)
  {
    Real tempD[C];
    Var<C> temp;
    for(int n=0; n<C;n++)
     {
       temp.m_ptr[n]=tempD+n;
       temp(n)=maskLoad(c(n), fixup); 
     }
    func(temp); // the end of the vector my have exception, but that's OK.
 
    for(int n=0; n<C; ++n)
    {
      maskStore(*(c.m_ptr[n]), *(temp.m_ptr[n]),fixup);
    }

  }
  // printf("size: %d, tripCount:%d fixup %d\n",size, tripCount, fixup);
}

template<typename Func, unsigned int C>
void forall(int size, Func func, Var<C> v1, Var<C> v2)
{
  //int tripCount=0;
  Var<C> c1 = v1, c2=v2;
  int i=0;
  for(; i<size-VLENGTH; i+=VLENGTH,++c1,++c2)
  {
    func(c1,c2); // tripCount++;
  }

  int fixup=size-VLENGTH*(size/VLENGTH);
  if(fixup>0)
  {
    Real tempD[2*C];
    Var<C> temp1, temp2;
    for(int n=0; n<C;n++)
     {
       temp1.m_ptr[n]=tempD+n;
       temp2.m_ptr[n]=tempD+C+n;
       temp1(n)=maskLoad(c1(n), fixup);
       temp2(n)=maskLoad(c2(n), fixup);
     }
    func(temp1, temp2); // the end of the vector my have exception, but that's OK.
    for(int n=0; n<C; ++n)
    {
      maskStore(*(c1.m_ptr[n]), *(temp1.m_ptr[n]),fixup);
      maskStore(*(c2.m_ptr[n]), *(temp2.m_ptr[n]),fixup);
    }
    /*
    std::cout<<"tempD "; for(int j=0; j<2*NUMCOMPS; ++j) std::cout<<tempD[j];
    std::cout<<"\n";
    std::cout<<"temp1 "<<temp1<<" c1 "<<c1<<std::endl;
    std::cout<<"temp2 "<<temp2<<" c2 "<<c2<<std::endl; 
    */
  }
  // printf("2 Vars. size: %d, tripCount:%d fixup %d\n",size, tripCount, fixup);
}
