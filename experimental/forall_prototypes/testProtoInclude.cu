
#define DIM 1
#define PROTO_CUDA 1
#include "../../include/Proto.H"
#include <cstdlib>
#include <cstdio>
#include <functional>
#include <iostream>

/* forall header material ============================ */
template<typename Func, typename... Rest>
__global__
void indexer(int begin, int end, Func body, Rest... a)
{
  int idx = threadIdx.x+blockIdx.x*blockDim.x;
  if(idx>=begin && idx<end)
  {
     body(idx, a...);
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

template<typename Func, typename... Rest>
inline
void
forall(int begin, int end, const Func& loop_body, Rest... a)
{
  constexpr int stride=8;
  const int blocks = (end-begin)/stride+1;
  indexer<<<stride, blocks>>>(begin, end, mapper(loop_body), a...);
}

//#define PROTO_KERNEL_START __device__ 
//#define PROTO_KERNEL_END(local_name, app_name) __device__ decltype(&local_name) app_name = local_name ;

//#define FORALL(begin, end, Func, ...)         \
//  const int block=(end-begin)/8+1; \
//  indexer<<<8, block>>>(begin, end, mapper(Func), __VA_ARGS__)

/*   End forall library code ============================ */

#define FORALL(begin, end, Func, ...)   \
  //__device__ decltype(&Func) dFunc = Func; \
  decltype(&Func) dFunc = Func; \
  static decltype(&Func) mFunc = mapper(dFunc);\
  indexer<<<8, (end-begin)/8+1>>>(begin, end, mFunc, __VA_ARGS__);

// User pointwise function
PROTO_KERNEL_START void initMultiF(int idx, int* a, int* b, int* c)
{
  a[idx]=0; b[idx]=idx; c[idx]=idx;
}
PROTO_KERNEL_END(initMultiF, initMulti)


// user application code

int main(int argc, char** argv) 
{
  constexpr int n = 16;

  int* dbuffer;  cudaMalloc(&dbuffer, 3*n*sizeof(int));
  int* aye=dbuffer, *bee=dbuffer+n, *cee=dbuffer+2*n;
  int hbuffer[3*n];
  int* a=hbuffer, *b=hbuffer+n, *c=hbuffer+2*n;

  forall(0, n, initMulti, aye, bee, cee);

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
 
 
  /* for my next trick.....
  forall(1,n-1,[] __device__ (int i, int* out, int* in){
      out[i]=in[i-1]-2*in[i]+in[i+1];},aye, bee);

  forall(0,n,[] __device__(int i, int* out, int* in){
      out[i]=in[i]*in[i];}, cee, bee); // asynchronous with previous function


  cudaMemcpy(hbuffer, dbuffer, 3*n*sizeof(int), cudaMemcpyDeviceToHost);
  cudaDeviceSynchronize();
  pass=true;
  for(int i=1; i<n-1; i++)
    {
      if(a[i]!=b[i-1]-2*b[i]+b[i+1] || c[i]!=b[i]*b[i]) pass=false;
    }

  if(pass) printf("PASS forall lambda\n");
  else     printf("FAIL forall lambda\n");

  */
  
  cudaFree(dbuffer);
  

  return 0;
}
