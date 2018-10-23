
#include <stdlib.h>
#include <stdio.h>
#include <functional>

//inline
__device__
void initMulti(int idx, int* a, int* b, int* c)
{
  a[idx]=0; b[idx]=idx; c[idx]=idx;
}

template<typename Func, typename... Rest>
__global__
void indexer(int begin, int end, Func body, Rest... a)
{
  int idx = threadIdx.x+blockIdx.x*blockDim.x;
  if(idx>=begin && idx<end)
  {
//does not work
    body(idx, a...);
//does work
//    initMulti(idx, a...);
  }
}

template<typename Func, typename... Rest>
void
forall(int begin, int end, Func loop_body, Rest... a)
{
  int stride=8;
  int blocks = (end-begin)/stride+1;
  indexer<<<stride, blocks>>>(begin, end, loop_body, a...);
}

__global__ 
void 
extractor_initMulti(void (**kernel)(int, (int*)...))
{

  *kernel = &initMulti;
}

int main(int argc, char** argv) 
{
  int n = 16;

  int* aye, *bee, *cee;
  cudaMalloc(&aye, n*sizeof(int));
  cudaMalloc(&bee, n*sizeof(int));
  cudaMalloc(&cee, n*sizeof(int));


  void (*h_ckernel1)(int* ...);
  void (**d_ckernel1)(int* ...);
  cudaMalloc(&d_ckernel1, sizeof(void *));
  extractor_initMulti<<<1,1>>>(d_ckernel1);
//  cudaMemcpy((void *)&h_ckernel1, (void *)d_ckernel1, sizeof(void *), cudaMemcpyDeviceToHost);

//  printf("made it to first forall\n");
//
//  //  forall<<< 1, 1>>> (0, n, &pointInitMultiple, aye, bee, cee);
//  forall(0, n, h_ckernel1, aye, bee, cee);
//
//  printf("going into cudaSynchronize \n");
//
//  //wait for gpu to finish before going back to cpu stuff
//  cudaDeviceSynchronize();
//
//  printf("out of cudaSynchronize \n");
//
//  int* a, *b, *c;
//  a = new int[n];
//  b = new int[n];
//  c = new int[n];
//  size_t bytes = n*sizeof(int);
//  cudaMemcpy(a, aye, bytes, cudaMemcpyDeviceToHost);
//  cudaMemcpy(b, bee, bytes, cudaMemcpyDeviceToHost);
//  cudaMemcpy(c, cee, bytes, cudaMemcpyDeviceToHost);
//  int a0 = a[0];
//  int b0 = b[0];
//  int c0 = c[0];
//
//
////  int a0 = aye[0];
////  int b0 = bee[0];
////  int c0 = cee[0];
//
//  printf(" a0= %i, b0= %i, c0 = %i\n", a0, b0, c0);
//
//
//  for(int i=0; i<n; ++i) 
//  {
//    //printf("i = %i, a= %i, b= %i, c = %i\n",i,  aye[i], bee[i], cee[i]);
//    printf("i = %i, a= %i, b= %i, c = %i\n",i,  a[i], b[i], c[i]);
//  }

  cudaFree(aye);
  cudaFree(bee);
  cudaFree(cee);
  cudaFree(d_ckernel1);
  return 0;
}
