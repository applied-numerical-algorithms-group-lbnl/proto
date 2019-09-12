#include <cstdio>

// nvcc -std=c++11 --expt-extended-lambda symbol.cu 

__device__ void testFunc(double& argA, int argB)
{
   argA*=argB;
}


//static function pointer or things break
__device__ decltype(&testFunc) testFuncF = testFunc;

__device__ void testFunc2(double& argA, int argB)
{
   argA+=argB;
}

__device__ decltype(&testFunc2) testFunc2F = testFunc2;


__global__ void kernel()
{
  printf("testFunc address from device: %p\n", &testFunc);
  printf("testFunc2 address from device: %p\n", &testFunc2);
}

template<typename Func>
__global__ void kernelArgs(Func func)
{
  double G=18;
  func(G , 12);
  printf(" function result=%f\n",G);
}
 
template<typename Func>
void forall(const Func& f)
{
  printf("function address from host: %p\n", &f);
  Func g(f);
  cudaMemcpyFromSymbol(&g, f, sizeof(g));
  printf("mapped function address from host: %p\n", g);
  kernelArgs<<<1,1>>>(g);  
}

template<typename Func>
void forlambda(const Func& f)
{
  kernelArgs<<<1,1>>>(f);
}

int main()
{
  printf("testFunc  address on host: %p\n",&testFunc);
  printf("testFunc2 address on host: %p\n",&testFunc2);
  fflush(stdout); 
  forall(testFuncF);
  cudaDeviceSynchronize();
  fflush(stdout); 

  forall(testFunc2F);
  cudaDeviceSynchronize();
  fflush(stdout); 

  forlambda([] __device__ (double& argA, int argB) { argA -= argB;});
  cudaDeviceSynchronize();
  fflush(stdout); 

  kernel<<<1,1>>>();
  cudaDeviceSynchronize();
  return 0;
}
