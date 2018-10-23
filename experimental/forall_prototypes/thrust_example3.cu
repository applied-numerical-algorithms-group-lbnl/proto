#include <thrust/transform_reduce.h>
#include <thrust/functional.h>
#include <thrust/device_vector.h>
#include <thrust/generate.h>
#include <thrust/reduce.h>
#include <thrust/functional.h>
#include <thrust/execution_policy.h>
#include <iostream>
template<typename T>
struct absolute_value 
{
  __host__ __device__ T operator()(const T &x) const
  {
    return x < T(0) ? -x : x;
  }
};


int main(void)
{
  int* d_data;
  size_t npts =6;
  cudaError err;
  cudaMalloc(&d_data, npts*sizeof(int));
  {
  thrust::device_ptr<int> devptr(d_data);

    
  err = cudaGetLastError();
  if (err != cudaSuccess)
  {
    fprintf(stderr, "cudaGetLastError() failed at %s:%i : %s\n",
            __FILE__, __LINE__, cudaGetErrorString(err));
  }
  
  int value = -42;
  thrust::fill(thrust::device, devptr, devptr + npts, value);
  err = cudaGetLastError();
  if (err != cudaSuccess)
  {
    fprintf(stderr, "cudaGetLastError() failed at %s:%i : %s\n",
            __FILE__, __LINE__, cudaGetErrorString(err));
  }
  }
/**/

  thrust::device_ptr<int> devptr(d_data);
  int result = thrust::transform_reduce(devptr, devptr + npts,
                                        absolute_value<int>(),
                                        0,
                                        thrust::maximum<int>());
  err = cudaGetLastError();
  if (err != cudaSuccess)
  {
    fprintf(stderr, "cudaGetLastError() failed at %s:%i : %s\n",
            __FILE__, __LINE__, cudaGetErrorString(err));
  }
  std::cout << "max value = " << result << std::endl;

/**/
  cudaFree(d_data);
  err = cudaGetLastError();
  if (err != cudaSuccess)
  {
    fprintf(stderr, "cudaGetLastError() failed at %s:%i : %s\n",
            __FILE__, __LINE__, cudaGetErrorString(err));
  }
  
}

