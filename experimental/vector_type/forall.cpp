
/*
 Testbed program for writing forall with gcc vector extensions

  compile with:

  >g++ -std=c++11 -mavx -O3 forall.cpp
  >clang++ -std=c++11 -mavx -O3 forall.cpp
 
or on Cori login node as
  >CC -std=c++11 -O3 forall.cpp

 ... moving on to testing for correctness on KNL node.

*/

#define NPOINTS 50

#include "PR_Real.h"
#include "forall.h"
#include "Euler.h"


int main(int argc, char* argv[])
{

    //  first just test to see if instrinsic functions can be called
    //  with Intel Vector types
   Real a(12.0*one), b(2.0*one), c(4.0*one);
   Real g=a*b+c;

   g=sqrt(g*M_PI);

   // build "BoxData" like things.
   constexpr int P=NPOINTS/VLENGTH+1;
   Real consFAB[P*NUMCOMPS];
   Real primFAB[P*NUMCOMPS];
   Real A[NUMCOMPS], B[NUMCOMPS];
   Var<NUMCOMPS> av, bv, cons, prim;
   for(int i=0; i<NUMCOMPS; i++)
   {
     cons.m_ptr[i]=consFAB+i*P;
     prim.m_ptr[i]=primFAB+i*P;
     av.m_ptr[i]=A+i;
     bv.m_ptr[i]=B+i;
   }

   //  Call some forall functions
   
   forall(NPOINTS, init, cons);

   forall(NPOINTS, consToPrim, prim, cons);


   //test results
   init(av);
   consToPrim(bv, av);
   std::cout<<"av: "<<av<<"\nbv: "<<bv<<std::endl;

   double* p = (double*)primFAB;
   bool pass=true;
   for(int n=0; n<NUMCOMPS; n++, p+=P*VLENGTH)
     for(int i=0; i<NPOINTS; i++)
       {
         if(p[i] != bv(n,0))
           {
             std::cout<<"prim["<<n<<","<<i<<"]="<<p[i]<<"\n";
             pass=false;
           }
       }
   if(pass) std::cout<<" vector init and consToPrim passed"<<std::endl;
   else std::cout<<" vector init and consToPrim failed"<<std::endl;

   for(int d=0; d<DIM; ++d)
     {
       dir=d; //get around variadic issue here.
       forall(NPOINTS, getFlux, prim, cons);
     }

   forall(NPOINTS, scaleInvDx, prim);
   
   return 0;
}
