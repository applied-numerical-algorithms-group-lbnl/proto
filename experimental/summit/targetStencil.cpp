
#include <cstdlib>
#include <cstdio>
#include <mpi.h>

#include "omp_stub.h"
#include <omp.h>
#include <vector>
#include <functional>
#include <type_traits>


#pragma omp declare target
void init(int idx, int* a, int* b, int* c)
{
  a[idx]=0;  b[idx]=idx ; c[idx]=idx;
}
#pragma omp end declare target
/*
template <typename Func, typename... T>
void forall(int begin, int end, Func loop_body, T&... t)
{
#pragma omp target is_device_ptr(t...)  // this line currently does not work. 
#pragma omp teams distribute parallel for simd
  for (int idx = begin; idx < end; ++idx) {
    loop_body(idx, t...);
  }
}
*/
template <typename Func, typename T>
void forall(int begin, int end, Func loop_body, T a, T b, T c)
{
  typename std::remove_reference<Func>::type body{loop_body};
#pragma omp target is_device_ptr(a,b,c) map(to:body)
#pragma omp teams distribute parallel for simd 
  for (int idx = begin; idx < end; ++idx) {
    body(idx, a,b,c);
  }
}
struct Stencil
{
  std::vector<int> m_coeff;
  std::vector<int>    m_offset;

  int* g_coeff=nullptr;
  int* g_offset=nullptr;
  std::size_t     n=0;
  int     g_devices=0;
  void close()
  {
    n = m_coeff.size();
    g_devices = omp_get_num_devices();
    int host =  omp_get_initial_device();
    if(g_devices >0)
      {
 
        for(int i=0; i<g_devices; ++i)
          {
            g_coeff = (int*)omp_target_alloc(n*sizeof(int),i);
            g_offset = (int*)omp_target_alloc(n*sizeof(int),i);
            omp_target_memcpy(g_coeff, m_coeff.data(),n*sizeof(int),
                              0,0,i, host);
            omp_target_memcpy(g_offset, m_offset.data(),n*sizeof(int),
                              0,0,i, host);
          }
      }
  }
  ~Stencil()
  {
    if(g_devices >0)
      {
        for(int i=0; i<g_devices; ++i)
          {
            omp_target_free(g_coeff, i);
            omp_target_free(g_offset,i);
          }
      }
  }
};

void apply(int begin, int end, const Stencil& stencil, int* src, int* dest)
{
  const int* coeff = (stencil.g_coeff==nullptr) ? stencil.m_coeff.data() :stencil.g_coeff;
  const int* offset= (stencil.g_coeff==nullptr) ? stencil.m_offset.data():stencil.g_offset;
  const int n = stencil.n;
#pragma omp target is_device_ptr(src, dest, coeff, offset)
#pragma omp parallel
#pragma teams distribute for simd 
  for (int idx = begin; idx < end; ++idx) {
    dest[idx]=0;
    for(int i=0; i<n; ++i){dest[idx]+=src[idx+offset[i]]*coeff[i];}
  }
}

int main(int argc, char* argv[]) {
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);
  int devices = omp_get_num_devices();
  int gpu = omp_get_default_device();
  int host =  omp_get_initial_device();
  int threads=1;
#pragma omp parallel
#pragma omp master
  {
    threads = omp_get_num_threads();
  }
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  printf("MPI rank %d has %d threads, %d devices. Host:%d Default device %d\n",
         rank, threads, devices, host, gpu);
  constexpr int n = 100;
  int* dbuffer = (int*)omp_target_alloc(3*n*sizeof(int),gpu);

  int *a=dbuffer, *b=dbuffer+n, *c=dbuffer+2*n;

  int hbuffer[3*n];
  int* ahost=hbuffer, *bhost=hbuffer+n, *chost=hbuffer+2*n;

  // using inline function
  //forall(0, n, init, a, b, c);
  forall(0,n,[](int i, int* a, int* b, int* c){a[i]=0; b[i]=c[i]=i;},
         a,b,c);
  printf("after init\n");
  // using lambda function
  forall(0, n,  [](int i, int* a, int* b, int* c){a[i]=b[i]*c[i];},
         a, b, c);

  printf("after lambda\n");
  omp_target_memcpy(hbuffer, dbuffer, sizeof(int)*3*n, 0,0,host,gpu);
  
  bool fail=false;
  for(int i=0; i<n; ++i)
    {
      if(ahost[i]!=i*i){ fail=true; printf("simple forall FAIL \n"); break;}
    }

  if(!fail) printf("lambda forall PASSED\n");

  //  Stencil testing

  Stencil simpleStencil;
  simpleStencil.m_coeff={1,-2,1};
  simpleStencil.m_offset={-1,0,1};
  simpleStencil.close();
  
  apply(1, n-1, simpleStencil, a, b);

  omp_target_memcpy(hbuffer, dbuffer, sizeof(int)*3*n, 0,0,host,gpu);
  fail = false;
  for(int i=1; i<n-1; ++i)
    {
      if(bhost[i]!=2){ fail=true; printf("simple 1D stencil FAIL \n"); break;}
    }
  
  if(!fail) printf("Stencil passed\n");

  omp_target_free(dbuffer,gpu);
  
  MPI_Finalize();
  return 0;
}
