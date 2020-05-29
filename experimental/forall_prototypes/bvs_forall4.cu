
#include <cstdlib>
#include <cstdio>
#include <functional>
#include <iostream>
#include <Proto_gpu.H>

/* forall header material ============================ */
template<typename Func, typename... Rest>
__global__
void indexer(int begin, int end, Func body, Rest... a)
{
  int idx = threadIdx.x+blockIdx.x*blockDim.x;
  if(idx<end)
  {
     body(idx, a...);
  }
}
// generic mapper to translate all function signatures
template<typename Func>
inline Func mapper(const Func& device_f)
{
  Func rtn(device_f); // trick needed for lambdas, since lambdas lack null constructors
  if (protoSuccess != protoMemcpyFromSymbol (&rtn, (const void*) device_f, sizeof (Func), 0, protoMemcpyDeviceToHost))
    printf ("FAILED to get SYMBOL\n");
  return rtn;
}

template<typename Func, typename... Rest>
inline
void
forall(int begin, int end, const Func& loop_body, Rest&&... a)
{
  constexpr int stride=8;
  const int blocks = (end-begin)/stride+1;
  protoLaunchKernel(indexer<Func,Rest...>, stride, blocks, begin, end, mapper(loop_body), std::forward<Rest>(a)...);
}

// User pointwise function
__device__ void initMultiF(int idx, int* a, int* b, int* c)
{
  a[idx]=0; b[idx]=idx; c[idx]=idx;
}
__device__ decltype(&initMultiF) initMulti = initMultiF;



// user application code

int main(int argc, char** argv) 
{
  constexpr int n = 16;

  int* dbuffer;  protoMalloc(&dbuffer, 3*n*sizeof(int));
  int* aye=dbuffer, *bee=dbuffer+n, *cee=dbuffer+2*n;
  int hbuffer[3*n];
  int* a=hbuffer, *b=hbuffer+n, *c=hbuffer+2*n;

  forall(0, n, initMulti, aye, bee, cee);

//  FORALL(0, n, initMultiF, aye, bee, cee);
  
  protoMemcpy(hbuffer, dbuffer, 3*n*sizeof(int), protoMemcpyDeviceToHost);

  bool pass=true;
  for(int i=0; i<n; ++i) 
    {
      if(a[i]!= 0 || b[i]!=i || c[i]!=i) pass=false;
      printf("i = %i, a= %i, b= %i, c = %i\n",i,  a[i], b[i], c[i]);
    }
  if(pass) printf("PASS init\n");
  else     printf("FAIL init\n");
 
 
  protoFree(dbuffer);
  

  return 0;
}
