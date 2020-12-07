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


template <class T>
class BoxDataMU
{
public:

  BoxDataMU(unsigned int a_size)
  {
    m_size = a_size;
    protoMalloc(m_deviceRawPtr, m_size*sizeof(T));
  }

  ~BoxDataMU()
  {
    protoFree(m_deviceRawPtr);
  }

  void setVal(T a_value)
  {
    T value = a_value;
    thrust::device_ptr<T> devptr = thrust::device_pointer_cast(m_deviceRawPtr);
    unsigned int nsize = m_size;

    thrust::fill(thrust::device, devptr, devptr+nsize,   value);

    protoError err = protoGetLastError();
    if (err != protoSuccess)
    {
      fprintf(stderr, "protoGetLastError() failed at %s:%i : %s\n",
              __FILE__, __LINE__, protoGetErrorString(err));
    }
  }

  T absMax()
  {
    thrust::device_ptr<T> devptr(m_deviceRawPtr);
    unsigned int nsize = m_size;
    T absmax1 = thrust::transform_reduce(devptr, devptr + nsize, 
                                         absolute_value<T>(),
                                         0,
                                         thrust::maximum<T>());

    protoError err = protoGetLastError();
    if (err != protoSuccess)
    {
      fprintf(stderr, "protoGetLastError() failed at %s:%i : %s\n",
              __FILE__, __LINE__, protoGetErrorString(err));
    }
    return absmax1;
  }

private:
  T* m_deviceRawPtr;
  unsigned int m_size;
  BoxDataMU();

};

int main(int argc, char **argv)
{

  size_t npts =8;
  BoxDataMU<int> data(npts);
  
  data.setVal(-42);

  int maxval = data.absMax();
  std::cout << "max abs value (should be 42) = " << maxval << std::endl;
  return 0;
}

