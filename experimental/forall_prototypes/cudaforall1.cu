
#include <stdlib.h>
#include <stdio.h>
#include <functional>


__device__
void pointInitMultiple(int idx,  int* a, int* b, int* c)
{
  a[idx]=0;
  b[idx]=idx;
  c[idx]=idx;
}


template <typename Func, typename... Rest >
__global__
void forall(int begin, int N, Func loop_body, Rest*... a)
{
  typename std::remove_reference<Func>::type body{loop_body};
  for(int idx = begin; idx < N; idx++)
    {
      body(idx, a ...);
    }
}

__global__
void pointInit2(int* a, int* b, int* c)
{
//  int index = threadIdx.x;
//  //number of threads in the block
//  int stride = blockDim.x;
  int idx = threadIdx.x;
  a[idx] = 0; 
  b[idx] = idx;
  c[idx] = idx;
}

template <typename Func, typename... Rest >
void forallbvs(int begin, int N, Func loop_body, Rest*... a)
{
  loop_body<<<1, N>>>(a...);
}

int main(int argc, char** argv) 
{
  int n = 2048;

  int* aye, *bee, *cee;
  cudaMalloc(&aye, n*sizeof(int));
  cudaMalloc(&bee, n*sizeof(int));
  cudaMalloc(&cee, n*sizeof(int));


  printf("made it to first forall\n");

  //  forall<<< 1, 1>>> (0, n, &pointInitMultiple, aye, bee, cee);
   forallbvs(0, n, &pointInit2, aye, bee, cee);

  printf("going into cudaSynchronize \n");

  //wait for gpu to finish before going back to cpu stuff
  cudaDeviceSynchronize();

  printf("out of cudaSynchronize \n");

int* a, *b, *c;
a = new int[n];
b = new int[n];
c = new int[n];
size_t bytes = n*sizeof(int);
  cudaMemcpy(a, aye, bytes, cudaMemcpyDeviceToHost);
  cudaMemcpy(b, bee, bytes, cudaMemcpyDeviceToHost);
  cudaMemcpy(c, cee, bytes, cudaMemcpyDeviceToHost);
  int a0 = a[0];
  int b0 = b[0];
  int c0 = c[0];


//  int a0 = aye[0];
//  int b0 = bee[0];
//  int c0 = cee[0];

  printf(" a0= %i, b0= %i, c0 = %i\n", a0, b0, c0);


  for(int i=0; i<n; ++i) 
    {
      //printf("i = %i, a= %i, b= %i, c = %i\n",i,  aye[i], bee[i], cee[i]);
      printf("i = %i, a= %i, b= %i, c = %i\n",i,  a[i], b[i], c[i]);
    }

  cudaFree(aye);
  cudaFree(bee);
  cudaFree(cee);

  return 0;
}
