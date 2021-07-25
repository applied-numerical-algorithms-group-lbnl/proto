#include <iostream>

struct test
{
  __device__ static 
  void gpu(double* ptr)
  {
    unsigned int idx = threadIdx.x + blockDim.x * blockIdx.x;
    ptr[idx] = 3;
  }
};

template<typename ker, typename... T>
__global__ 
void launcher(T... args)
{
  ker::gpu(args...);
}

int main()
{
  unsigned int size = 4;
  double* myptr;
  cudaMalloc(&myptr,size*sizeof(double));
  launcher<test,double*><<<1,size>>>(myptr);

  double *host = new double[size];
  cudaMemcpy(host,myptr,size*sizeof(double),cudaMemcpyDeviceToHost);

  for(int i = 0; i<size ;i++)
    std::cout << host[i] << std::endl;

  return 0;
}
