#include <thrust/transform_reduce.h>
#include <thrust/functional.h>
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
  int data[6] = {-1, 0, -2, -2, 1, -3};
  int result = thrust::transform_reduce(data, data + 6,
                                        absolute_value<int>(),
                                        0,
                                        thrust::maximum<int>());
  std::cout << "max value = " << result << std::endl;
  
}

