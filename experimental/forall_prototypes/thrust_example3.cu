#include <thrust/transform_reduce.h>
#include <thrust/functional.h>
#include <thrust/device_vector.h>
#include <thrust/generate.h>
#include <thrust/reduce.h>
#include <thrust/functional.h>
#include <thrust/execution_policy.h>
#include <iostream>
#include <Proto_gpu.H>
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
  protoError err;
  protoMalloc(d_data, npts*sizeof(int));
  {
  thrust::device_ptr<int> devptr(d_data);

    
  err = protoGetLastError();
  if (err != protoSuccess)
  {
    fprintf(stderr, "protoGetLastError() failed at %s:%i : %s\n",
            __FILE__, __LINE__, protoGetErrorString(err));
  }
  
  int value = -42;
  thrust::fill(thrust::device, devptr, devptr + npts, value);
  err = protoGetLastError();
  if (err != protoSuccess)
  {
    fprintf(stderr, "protoGetLastError() failed at %s:%i : %s\n",
            __FILE__, __LINE__, protoGetErrorString(err));
  }
  }
/**/

  thrust::device_ptr<int> devptr(d_data);
  int result = thrust::transform_reduce(devptr, devptr + npts,
                                        absolute_value<int>(),
                                        0,
                                        thrust::maximum<int>());
  err = protoGetLastError();
  if (err != protoSuccess)
  {
    fprintf(stderr, "protoGetLastError() failed at %s:%i : %s\n",
            __FILE__, __LINE__, protoGetErrorString(err));
  }
  std::cout << "max value = " << result << std::endl;

/**/
  protoFree(d_data);
  err = protoGetLastError();
  if (err != protoSuccess)
  {
    fprintf(stderr, "protoGetLastError() failed at %s:%i : %s\n",
            __FILE__, __LINE__, protoGetErrorString(err));
  }
  
}

