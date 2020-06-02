#include <cstdio>
#include <cstdint>
#include <Proto_gpu.H>

__device__ int f() { return 42 ; }

__device__ decltype(&f) f_ptr = f ;

template<typename Func>
__global__ void kernel ( Func func )
{
    int k = func() ;
    printf ("%d\n", k) ;
}

// generic mapper to translate all function signatures
template<typename Func>
inline Func mapper(const Func& device_f)
{
  Func rtn(device_f);
  if (protoSuccess != protoMemcpyFromSymbol (&rtn, (const void *)device_f, sizeof (Func), 0, protoMemcpyDeviceToHost))
    printf ("FAILED to get SYMBOL\n");
  return rtn;
}
  
int main ()
{
 
  protoLaunchKernel(kernel, 1, 1, mapper(f_ptr));
  protoDeviceSynchronize();
}
