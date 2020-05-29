
#include <stdlib.h>
#include <stdio.h>
#include <functional>
#include "../../include/Proto_gpu.H"

__global__
void init(int n, int* a, int* b, int* c, int* d)
{
  //index of current thread
  int index = threadIdx.x;
  //number of threads in the block
  int stride = blockDim.x;
  for(int idx = index; idx < n; idx += stride)
    {
      a[idx]=0;
      b[idx]=idx;
      c[idx]=idx;
      d[idx]=0;
    }

}

__global__
void pointInit(int idx,  int* d)
{
  d[idx]=0;
}

__global__
void pointInitMultiple(int idx,  int* a, int* b, int* c, int* d )
{
  a[idx]=0;
  b[idx]=idx;
  c[idx]=idx;
  d[idx]=0;
}

__global__
void incrementAwithBplusC(int idx,  int* a, int* b, int* c)
{
  a[idx] += b[idx] + c[idx];
}


template <typename Func, typename... Rest >
__global__
void forall(int begin, int end, Func loop_body, Rest*... a)
{
  //index of current thread
  int index = threadIdx.x;
  //number of threads in the block
  int stride = blockDim.x;

  //typename std::remove_reference<Func>::type body{loop_body};
  int distance = end-begin;
  int start = begin + index;
  int n = distance;
  for(int idx = start; idx < n; idx += stride)
    {
      loop_body(idx, a ...);
    }
}

template <typename Func, typename... Rest >
void trickLaunchKernelForall(int n, int m, int begin, int end, Func loop_body, Rest*... a)
{
	protoLaunchKernel(forall<Func,Rest...>, n, m, begin, end, loop_body, a...);
}

int main(int argc, char** argv) 
{
  int n = 2048;
  /**
     CUDA GPUs have many parallel processors grouped into Streaming
     Multiprocessors, or SMs. Each SM can run multiple concurrent
     thread blocks. As an example, a Tesla P100 GPU based on the
     Pascal GPU Architecture has 56 SMs, each capable of supporting up
     to 2048 active threads. To take full advantage of all these
     threads, I should launch the kernel with multiple thread blocks.

     By now you may have guessed that the first parameter of the execution
     configuration specifies the number of thread blocks. Together, the
     blocks of parallel threads make up what is known as the grid. Since I
     have N elements to process, and 256 threads per block, I just need to
     calculate the number of blocks to get at least N threads. I simply
     divide N by the block size (being careful to round up in case N is not
     a multiple of blockSize).
   **/
  int blockSize = 256;
  int numBlocks = (n + blockSize-1)/blockSize;

  int* aye, *bee, *cee, *dee;
  protoMallocManaged(&aye, n*sizeof(int));
  protoMallocManaged(&bee, n*sizeof(int));
  protoMallocManaged(&cee, n*sizeof(int));
  protoMallocManaged(&dee, n*sizeof(int));

  protoLaunchKernel(init, numBlocks, blockSize, n, aye, bee, cee, dee);


  printf("made it to first forall\n");
  trickLaunchKernelForall(numBlocks, blockSize, 0, n, &pointInit, dee);

  printf("made it to second forall\n");
  trickLaunchKernelForall(numBlocks, blockSize, 0, n, &pointInitMultiple, aye, bee, cee, dee);

  printf("made it to third forall\n");
  trickLaunchKernelForall(numBlocks, blockSize, 0, n, &incrementAwithBplusC, aye, bee, cee);

  printf("going into cudaSynchronize \n");
  //wait for gpu to finish before going back to cpu stuff
  protoDeviceSynchronize();

  printf("out of cudaSynchronize \n");
  int* a, *b, *c;
  a = new int[n];
  b = new int[n];
  c = new int[n];
  size_t bytes = n*sizeof(int);
  protoMemcpy(a, aye, bytes, protoMemcpyDeviceToHost);
  protoMemcpy(b, bee, bytes, protoMemcpyDeviceToHost);
  protoMemcpy(c, cee, bytes, protoMemcpyDeviceToHost);

  bool pass = true;
  for(int i=0; i<n; ++i) 
    {
      if (a[i] != 2*i) pass = false;
      if(!pass)
	{
	  printf("a= %i, b= %i, c = %i\n", a[i], b[i], c[i]);
	}
    }
  if(pass)
    {
      printf("PASS\n");
    }
  else
    {
      printf("FAIL\n");
    }

  protoFree(aye);
  protoFree(bee);
  protoFree(cee);
  protoFree(dee);

  delete[] a;
  delete[] b;
  delete[] c;

  return 0;
}
