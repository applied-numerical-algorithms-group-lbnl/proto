
/*
 Testbed program for writing forall with gcc vector extensions

  compile with:

  >g++ -std=c++11 -mavx -O3 forall.cpp
  >clang++ -std=c++11 -mavx -O3 forall.cpp
 
or on Cori login node as
  >CC -std=c++11 -O3 forall.cpp

 ... moving on to testing for correctness on KNL node.

*/

#define NPOINTS (32*32*32)
#define NBOXES 50
#define ITERATIONS 25

#include "PR_Real.h"
#include "forall.h"
#include "Euler.h"
#include <vector>
#include <ctime>
#include <ratio>
#include <chrono>

 using namespace std::chrono;

struct BData
{
  Real m_data[(NPOINTS/VLENGTH+1)*NUMCOMPS];
  Var<NUMCOMPS> start()
  {
    Var<NUMCOMPS> rtn;
    for(int i=0; i<NUMCOMPS; i++)
      {
        rtn.m_ptr[i]=m_data+i*(NPOINTS/VLENGTH+1);
      }
    return rtn;
  }
};

int main(int argc, char* argv[])
{
  

   // build "BoxData" like things.
   constexpr int P=NPOINTS/VLENGTH+1;
   Real consFAB[P*NUMCOMPS];
   Real primFAB[P*NUMCOMPS];
   Real A[NUMCOMPS], B[NUMCOMPS];
   Var<NUMCOMPS> cons, prim;


   std::vector<BData> consV(NBOXES);
   std::vector<BData> primV(NBOXES);

   for(int i=0; i<NBOXES; i++)
     {
       auto cons=consV[i].start();
       auto prim=primV[i].start();
       forall(NPOINTS, init, cons);
     }
   high_resolution_clock::time_point t1 = high_resolution_clock::now();
   //  Call some forall functions
   for(int t=0; t<ITERATIONS; t++)
     {

       for(int i=0; i<NBOXES; i++)
         {
           auto cons=consV[i].start();
           auto prim=primV[i].start();
           
           forall(NPOINTS, consToPrim, prim, cons);
           
           for(int d=0; d<DIM; ++d)
             {
               dir=d; //get around variadic issue here.
               forall(NPOINTS, getFlux, prim, cons);
             }
           
           forall(NPOINTS, scaleInvDx, prim);
         }
     }
   high_resolution_clock::time_point t2 = high_resolution_clock::now();

   duration<double> time_span = duration_cast<duration<double>>(t2 - t1);

   std::cout << "It took me " << time_span.count() << " seconds.";
   std::cout << std::endl;
   int flops = ITERATIONS*NPOINTS*NBOXES*(7+DIM*10+1);
   std::cout<<"flops: "<<flops<<" GFLOPS/s: "<<(flops/time_span.count())/1000000000;
   std::cout<<std::endl;
   return 0;
}
