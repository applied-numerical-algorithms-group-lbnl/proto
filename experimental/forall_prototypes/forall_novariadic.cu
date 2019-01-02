
#include <cstdlib>
#include <cstdio>
#include <functional>
#include <iostream>

/* forall header material ============================ */
template<typename Func>
__global__
void indexer(int begin, int end, Func body, int* a, int* b, int* c)
{
  int idx = threadIdx.x+blockIdx.x*blockDim.x;
  if(idx<end)
  {
    body(idx, a, b, c);
  }
}
// generic mapper to translate all function signatures
template<typename Func>
inline Func mapper(const Func& device_f)
{
  Func rtn(device_f); // trick needed for lambdas, since lambdas lack null constructors
  if (cudaSuccess != cudaMemcpyFromSymbol (&rtn, device_f, sizeof (Func)))
    printf ("FAILED to get SYMBOL\n");
  return rtn;
}

template<typename Func>
inline
void
forall(int begin, int end, int* a, int* b, int* c)
{
  constexpr int stride=8;
  const int blocks = (end-begin)/stride+1;
  indexer<<<stride, blocks>>>(begin, end, &Func::op, a, b, c);
}

// User pointwise function
__device__ void initMultiF(int idx, int* a, int* b, int* c)
{
  a[idx]=0; b[idx]=idx; c[idx]=idx;
}
struct initMulti { static __device__ void op(int idx, int* a, int* b, int* c)
  { return initMultiF(idx, a, b, c);}};


// user application code

int main(int argc, char** argv) 
{
  constexpr int n = 16;

  int* dbuffer;  cudaMalloc(&dbuffer, 3*n*sizeof(int));
  int* aye=dbuffer, *bee=dbuffer+n, *cee=dbuffer+2*n;
  int hbuffer[3*n];
  int* a=hbuffer, *b=hbuffer+n, *c=hbuffer+2*n;

  forall<initMulti>(0, n, aye, bee, cee);

//  FORALL(0, n, initMultiF, aye, bee, cee);
  
  cudaMemcpy(hbuffer, dbuffer, 3*n*sizeof(int), cudaMemcpyDeviceToHost);

  bool pass=true;
  for(int i=0; i<n; ++i) 
    {
      if(a[i]!= 0 || b[i]!=i || c[i]!=i) pass=false;
      printf("i = %i, a= %i, b= %i, c = %i\n",i,  a[i], b[i], c[i]);
    }
  if(pass) printf("PASS init\n");
  else     printf("FAIL init\n");
 
 
  cudaFree(dbuffer);
  

  return 0;
}
