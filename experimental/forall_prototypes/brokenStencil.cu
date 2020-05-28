
#include <cstdlib>
#include <cstdio>
#include <functional>
#include <iostream>
#include <vector>

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
  if (protoSuccess != protoMemcpyFromSymbol (&rtn, device_f, sizeof (Func)))
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

#define PROTO_KERNEL_START __device__ 
#define PROTO_KERNEL_END(local_name, app_name) __device__ decltype(&local_name) app_name = local_name ;


// User pointwise function
PROTO_KERNEL_START void initMultiF(int idx, int* a, int* b, int* c)
{
  a[idx]=idx*idx; b[idx]=0; c[idx]=0;
}
PROTO_KERNEL_END(initMultiF, initMulti)

struct Stencil
{
  //lives on host
  std::vector<int> m_coeff;
  std::vector<int>    m_offset;

  //lives on device
  int* g_coeff=nullptr;
  int* g_offset=nullptr;
  std::size_t     n=0;

  void close()
  {
    n = m_coeff.size();
    size_t memsize = n*sizeof(int);
 
    protoMalloc(&g_coeff,  memsize);
    protoMalloc(&g_offset, memsize);
    protoMemcpy(g_coeff ,  m_coeff.data(), memsize, protoMemcpyHostToDevice);
    protoMemcpy(g_offset, m_offset.data(), memsize, protoMemcpyHostToDevice);
  }
  ~Stencil()
  {
    protoFree(g_coeff );
    protoFree(g_offset);
  }
};

PROTO_KERNEL_START void applyInnerLoopF(int idx, int* nptr, int* src, int* dst, int* coeff, int* offset)
{
  int n = *nptr;
  for(int i=0; i<n; ++i)
  {
    dst[idx]+=src[idx+offset[i]]*coeff[i];
  }
}
PROTO_KERNEL_END(applyInnerLoopF, applyInnerLoop)

PROTO_KERNEL_START void addFourInnerLoopF(int idx, int* nptr, int* src, int* dst, int* coeff, int* offset)
{
  int n = *nptr;
  for(int i=0; i<n; ++i)
  {
    dst[idx]+= 4;
  }
}
PROTO_KERNEL_END(addFourInnerLoopF, addFourInnerLoop)

void apply(int begin, int end, Stencil& stencil, int* src, int* dest)
{
  int* coeff = (stencil.g_coeff==nullptr) ? stencil.m_coeff.data() :stencil.g_coeff;
  int* offset= (stencil.g_coeff==nullptr) ? stencil.m_offset.data():stencil.g_offset;
  int n = stencil.n;
//  constexpr int stride=8;
//  const int blocks = (end-begin)/stride+1;
  for (int idx = begin; idx < end; ++idx) 
  {
    forall(0, n, addFourInnerLoop, idx, &n, src, dest, coeff, offset);
  }
//  indexer<<<stride, blocks>>>(begin, end, mapper(applyInnerLoop), &n, src, dest, coeff, offset);
//  indexer<<<stride, blocks>>>(begin, end, mapper(addFourInnerLoop), &n, src, dest, coeff, offset);
}

int main(int argc, char* argv[]) 
{
  constexpr int n = 16;
  int* devbuffer;

  protoMalloc(&devbuffer, 3*n*sizeof(int));
  int hostbuffer[3*n];
  //bvs-- evil genius at work
  int *adev=devbuffer, *bdev=devbuffer+n, *cdev=devbuffer+2*n;
  int* ahost=hostbuffer; 
  int* bhost=hostbuffer + n;
  int* chost=hostbuffer+2*n; 

  forall(0, n, initMulti, adev, bdev, cdev);
  protoDeviceSynchronize();
  printf("after init\n");
  protoMemcpy(hostbuffer, devbuffer, sizeof(int)*3*n, protoMemcpyDeviceToHost);
  printf("after cudamemcopy \n");

  for(int i=1; i<n-1; ++i)
    {
      printf("i = %d, ahost[i]=%d, bhost[i]=%d, chost[i] = %d \n", i, ahost[i],bhost[i], chost[i]);
    }

  Stencil simpleStencil;
  simpleStencil.m_coeff={1,-2,1};
  simpleStencil.m_offset={-1,0,1};
  simpleStencil.close();

  printf("after stencil closure \n");

  apply(1, n-1, simpleStencil, adev, bdev);
   
  protoDeviceSynchronize();
  printf("after stencil apply\n");
  protoMemcpy(hostbuffer, devbuffer, sizeof(int)*3*n, protoMemcpyDeviceToHost);

  printf("after cudamemcopy \n");
  for(int i=1; i<n-1; ++i)
    {
      printf("i = %d, ahost[i]=%d, bhost[i]=%d, chost[i] = %d \n", i, ahost[i],bhost[i], chost[i]);
    }

  bool fail = false;
  for(int i=1; i<n-1; ++i)
    {
      int value = bhost[i];
      printf("i = %d, bhost[i] = %d \n", i, value);
      if(bhost[i]!=2){ fail=true; printf("simple 1D stencil FAIL \n"); break;}
    }
  
  if(!fail) printf("Stencil passed\n");

  protoFree(devbuffer);
  
  return 0;
}
