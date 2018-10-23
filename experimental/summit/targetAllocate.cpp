
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include "omp_stub.h"
#include <omp.h>

template <typename Func>
void forall(int begin, int end, Func loop_body)
{
  #pragma teams distribute parallel for simd
  for (int idx = begin; idx < end; ++idx) {
    loop_body(idx);
  }
}

int main(int argc, char* argv[]) {
  int provided;
  printf("Starting program\n");
//  MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);
  MPI_Init(&argc, &argv);
  int devices = omp_get_num_devices();
  int gpu = omp_get_default_device();
  int host =  omp_get_initial_device();
  int threads = omp_get_num_threads();
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  printf("MPI rank %d has %d threads and %d devices: host=%d gpu=%d\n",rank, threads, devices,host, gpu);

  constexpr int n = 10;
  int* dbuffer = (int*)omp_target_alloc(3*n*sizeof(int),gpu);

  int *a=dbuffer, *b=dbuffer+n, *c=dbuffer+2*n;

  int hbuffer[3*n];
  int* ahost=hbuffer, *bhost=hbuffer+n, *chost=hbuffer+2*n;


 printf("set val loop\n"); 
#pragma omp target is_device_ptr(a,b,c) device(gpu) 
#pragma omp teams distribute parallel for simd
  for(int i=0; i<n; ++i) {
    a[i] = 0; b[i] = c[i] = i;
  }

  printf("set val forall lambda\n");
#pragma omp target is_device_ptr(a,b,c) device(gpu)
  forall(0,n,[=](int i) { a[i]=0; b[i]=c[i]=i;});

#pragma omp target is_device_ptr(a,b,c) device(gpu)
  forall(0, n, [&] (int i) {
    a[i] += b[i] + c[i];
  });

  printf("omp_target_memcpy\n");
  omp_target_memcpy(hbuffer, dbuffer, sizeof(int)*3*n, 0,0,host,gpu);
  
  if(rank==0)
   {
     printf("\n");
     for(int i=0; i<n; ++i) printf("%d ",ahost[i]);
     printf("\n");
   }
     
  for(int i=0; i<n; ++i) 
    if (ahost[i] != 2*i) {printf("%d FAIL\n",rank); break;}

  omp_target_free(dbuffer,gpu);
  
  MPI_Finalize();
  return 0;
}
