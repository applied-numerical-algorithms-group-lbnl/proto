
#include <stdlib.h>
#include <functional>
#include "omp_stub.h"
#include <omp.h>

#pragma omp declare target 
static void init(int idx, int* a)
{
  a[idx]=0;
}

#pragma omp end declare target

template <typename Func>
void forall(int begin, int end, Func loop_body, int* a)
{
  //typename std::remove_reference<Func>::type body{loop_body};
  int distance = end-begin;
#pragma omp target teams distribute parallel for is_device_ptr(a,loop_body) 
  for (int idx = 0; idx < distance; ++idx) {
    loop_body(idx, a);
  }
}
template <typename Func>
void forall(int begin, int end, Func loop_body)
{
  typename std::remove_reference<Func>::type body{loop_body};
  int distance = end-begin;
  #pragma omp target teams distribute parallel for map(to:body)
  for (int idx = 0; idx < distance; ++idx) {
    body(idx);
  }
}
int main() {
  const int n = 10;
  int a[n], b[n], c[n];
  for(int i=0; i<n; ++i) {
    a[i] = 0; b[i] = c[i] = i;
  }
  int gpu = omp_get_default_device();
  int* a_d = (int*)omp_target_alloc(n*sizeof(int), gpu);
#pragma omp target teams distribute parallel for is_device_ptr(a_d)
  for(int i=0; i<n; ++i) { init(i, a_d);}

  printf("made it to first forall\n");
  
  forall(0, n, &init, a_d); // compiles but throws runtime error

  printf("made it to second forall\n");
  #pragma omp target enter data map(to: a, b[:n], c)
  forall(0, n, [&] (int i) {
    a[i] += b[i] + c[i];
  });
  #pragma omp target exit data map(from: a[:n]) map(release:b[:n], c[:n])
 
  for(int i=0; i<n; ++i) 
    if (a[i] != 2*i) printf("FAIL\n");
  printf("PASS\n");
  return 0;
}
