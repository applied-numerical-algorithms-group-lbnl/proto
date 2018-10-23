
#include <stdlib.h>
#include <stdio.h>
#include <functional>


void pointInitMultiple(int idx,  int* a, int* b, int* c)
{
  a[idx]=0;
  b[idx]=idx;
  c[idx]=idx;
}

void incrementAwithBplusC(int idx,  int* a, int* b, int* c)
{
  a[idx] += b[idx] + c[idx];
}


template <typename Func, typename... Rest >
void forall(int begin, int N, Func loop_body, Rest*... a)
{
//  //index of current thread
//  int index = threadIdx.x;
//  //number of threads in the block
//  int stride = blockDim.x;
//
//  //typename std::remove_reference<Func>::type body{loop_body};
//  //  int distance = end-begin;
//  int start = begin + index;
//  //  int n = distance;
//  //  for(int idx = start; idx < end; idx += stride)
  for(int idx = begin; idx < N; idx++)
    {
      loop_body(idx, a ...);
    }
}

int main(int argc, char** argv) 
{
  int n = 16;
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

  numBlocks = 1;
  blockSize = 1;

  int *aye = new int[n];
  int *bee = new int[n];
  int *cee = new int[n];

  printf("made it to first forall\n");
  forall(0, n, &pointInitMultiple, aye, bee, cee);

  for(int i=0; i<n; ++i) 
    {
	  printf("i = %i, a= %i, b= %i, c = %i\n",i,  aye[i], bee[i], cee[i]);
    }

  delete[](aye);
  delete[](bee);
  delete[](cee);

  return 0;
}
