
#include <cstdlib>
#include <cstdio>
#include <functional>
#include <iostream>
#include <vector>
/* stencil apply  header material ============================ */
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
    printf("stencil.size() = %d\n", n);
    for(int isten = 0; isten < n; isten++)
    {
      printf("coeff[%d])  = %d\n",isten, m_coeff[isten]);
      printf("offset[%d]) = %d\n",isten, m_offset[isten]);
    }
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

__global__
void stencilIndexer(int a_begin, int a_end, int a_n, int* a_src, int* a_dst,int* coeff,int* offset)
{
  int idx = a_begin + threadIdx.x+blockIdx.x*blockDim.x;
  int n = a_n;
  printf("inside indexer n = %d, begin = %d, end = %d, idx = %d\n", n, a_begin, a_end, idx);
  if((idx >= a_begin) && (idx< a_end))
  {
    for(int isten = 0; isten < n; isten++)
    {
      a_dst[idx] += a_src[idx+offset[isten]]*coeff[isten];
    }
  }
}

__global__
void addFourIndexer(int a_begin, int a_end, int a_n, int* a_src, int* a_dst)
{
  int idx = a_begin + threadIdx.x+blockIdx.x*blockDim.x;
  int n = a_n;
  printf("inside indexer n = %d, begin = %d, end = %d, idx = %d\n", n, a_begin, a_end, idx);
  if((idx >= a_begin) && (idx< a_end))
  {
    for(int isten = 0; isten < n; isten++)
    {
      a_dst[idx] += 4;
    }
  }
}



void
apply(int begin, int end, Stencil& a_stencil, int* a_src, int* a_dst)
{
//  constexpr int stride=8;
//  const int blocks = (end-begin)/stride+1;
  int n = a_stencil.n;
  printf("going into indexer with n=  %d, begin = %d, end = %d\n", n, begin, end);
  constexpr int stride=8;
  int blocks = (end-begin)/stride+1;
  int* coeff  =  a_stencil.g_coeff;
  int* offset =  a_stencil.g_offset;
//  addFourIndexer<<<stride, blocks>>>(begin, end, n, a_src, a_dst);
  stencilIndexer<<<stride, blocks>>>(begin, end, n, a_src, a_dst, coeff, offset);
}


//void hostAddFourIndexer(int idx, int a_begin, int a_end, Stencil& a_stencil, int* a_src, int* a_dst)
//{
//  int n = a_stencil.n;
//  printf("n = %d, begin = %d, end = %d, idx = %d\n", n, a_begin, a_end, idx);
//  if((idx >= a_begin) && (idx< a_end))
//  {
//    for(int isten = 0; isten < a_stencil.n; isten++)
//    {
//      a_dst[idx] += 4;
//    }
//  }
//}
//
//
//void
//hostApply(int begin, int end, Stencil& a_stencil, int* a_src, int* a_dst)
//{
//  int n = a_stencil.n;
//  printf("host going into indexer with n=  %d, begin = %d, end = %d\n", n, begin, end);
//  for(int idx = begin; idx <= end; idx++)
//  {
//    hostAddFourIndexer(idx, begin, end, a_stencil, a_src, a_dst);
//  }
//}

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
Func mapper(const Func& device_f)
{
  Func rtn(device_f); // trick needed for lambdas, since lambdas lack null constructors
  protoError_t err = protoMemcpyFromSymbol(&rtn, device_f, sizeof(Func), 0, protoMemcpyDeviceToHost);
  if (err != protoSuccess)
  {
    printf("FAILED to get SYMBOL\n");
    fprintf(stderr, "protoMemcpyFromSymbol() failed at %s:%i : %s\n",
            __FILE__, __LINE__, protoGetErrorString(err));
  }
  return rtn;
}

template<typename Func, typename... Rest>

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


//  hostApply(1, n-1, simpleStencil, ahost, bhost);
//  printf("after host apply  \n");
//  for(int i=1; i<n-1; ++i)
//    {
//      printf("i = %d, ahost[i]=%d, bhost[i]=%d, chost[i] = %d \n", i, ahost[i],bhost[i], chost[i]);
//    }

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
