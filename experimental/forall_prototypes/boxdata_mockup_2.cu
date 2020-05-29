#define DIM 2
#define CUDA_DECORATION __host__ __device__

#include <thrust/transform_reduce.h>
#include <thrust/functional.h>
#include <thrust/device_vector.h>
#include <thrust/generate.h>
#include <thrust/reduce.h>
#include <thrust/functional.h>
#include <thrust/execution_policy.h>
#include <iostream>
#include <Proto_gpu.H>

class PointMU
{
public:
  CUDA_DECORATION inline PointMU(const int a_tuple[DIM])
  {
    for (int j = 0; j < DIM; j++)
    {
      m_tuple[j]=a_tuple[j];
    }
  }

  inline static PointMU Zeros()
  {
    int zeros[DIM];
    for (int k = 0 ; k < DIM; k++)
    {
      zeros[k] = 0;
    }
    return PointMU(zeros);
  }

  inline static PointMU Ones(int a_scale=1)
  {
    int ones[DIM];
    for (int k = 0 ; k < DIM; k++)
    {
      ones[k] = a_scale;
    }
    return PointMU(ones);
  } 

  int operator[](unsigned int a_dir) const
  {
    return m_tuple[a_dir];
  }

  PointMU()
  {
  }
private:


  int m_tuple[DIM]; 
};

class BxMU
{
public:
  BxMU(const PointMU& a_low, const PointMU& a_high)
  {
    m_low  = a_low;
    m_high = a_high;
  }

  inline static BxMU Cube(int a_size)
  {
    return BxMU(PointMU::Zeros(),PointMU::Ones(a_size-1));
  }
  inline size_t size() const
  {
    size_t retval = m_high[0]-m_low[0] + 1;
    for(int idir = 1; idir < DIM; idir++)
    {
      size_t dimsize = m_high[idir]-m_low[idir] + 1;
      retval *= dimsize;
    }
    return retval;
  }
  BxMU()
  {
    m_low  = PointMU::Ones();
    m_high = PointMU::Zeros();
  }
  
private:
  PointMU m_low; ///< Point object containing the lower bounds of the Bx.
  PointMU m_high;  ///< Point object containing the upper bounds of the Bx.
};

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

  BoxDataMU(const BxMU& a_box)
  {
    m_box  = a_box;
    m_size = a_box.size();
    protoMalloc(&m_deviceRawPtr, m_size*sizeof(T));
    protoError err = protoGetLastError();
    if (err != protoSuccess)
    {
      fprintf(stderr, "protoGetLastError() failed at %s:%i : %s\n",
              __FILE__, __LINE__, protoGetErrorString(err));
    }
    m_deviceData = ::std::shared_ptr<T>(m_deviceRawPtr, [](T* p){protoFree(p);});
  }

  ~BoxDataMU()
  {
    protoFree(m_deviceRawPtr);
  }

  void setVal(T a_value)
  {
    T value = a_value;
    thrust::device_ptr<T> devptr = thrust::device_pointer_cast(m_deviceRawPtr);
    size_t nsize = m_size;

    thrust::fill(thrust::device, devptr, devptr+nsize,   value);

    protoError err = protoGetLastError();
    if (err != protoSuccess)
    {
      fprintf(stderr, "protoGetLastError() failed at %s:%i : %s\n",
              __FILE__, __LINE__, protoGetErrorString(err));
    }
  }

  T absMax() const
  {
    thrust::device_ptr<T> devptr(m_deviceRawPtr);
    size_t nsize = m_size;
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
  
  inline std::size_t size() const {return m_box.size();};

private:
  T* m_deviceRawPtr;
  BxMU m_box;
  ::std::shared_ptr<T>    m_deviceData; ///< Data array
  size_t m_size;
  BoxDataMU();

};

int main(int argc, char **argv)
{

  size_t npts =8;
  BxMU domain = BxMU::Cube(npts);
  BoxDataMU<int> data(domain);
  
  data.setVal(-42);

  int maxval = data.absMax();
  std::cout << "max abs value (should be 42) = " << maxval << std::endl;
  return 0;
}

