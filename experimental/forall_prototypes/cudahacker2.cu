#include <stdio.h>

template <typename T>
__global__ void ckernel1(T *data){

  int my_val = (int)(*data+1);
  printf("hello: %d \n", my_val);
}
template <typename TFunc, typename... TArgs>
__global__ void Test(TFunc func, int count, TArgs... args)
{
#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ >= 350)
  (*func)<< <1, 1 >> >(args...);
#else
  printf("What are you doing here!?\n");
#endif
}

template <typename... TArgs>
__host__ void Iterate(void(*kernel)(TArgs...), const int sysInfo, int count, TArgs... args)
{
  if(sysInfo >= 350)
  {
    printf("Iterate on GPU\n");
    Test << <1, 1 >> >(kernel, count, args...);
  }
  else
  {
    printf("Iterate on CPU\n");
    Test << <1, 1 >> >(kernel, count, args...);
  }
}

template <typename T>
__global__ void extractor(void (**kernel)(T *)){

  *kernel = ckernel1<T>;
}

template <typename T>
void run_test(T init)
{

  void (*h_ckernel1)(T *);
  void (**d_ckernel1)(T *);
  T *d_data;
  cudaMalloc(&d_ckernel1, sizeof(void *));
  cudaMalloc(&d_data, sizeof(T));
  cudaMemcpy(d_data, &init, sizeof(T), cudaMemcpyHostToDevice);
  extractor<<<1,1>>>(d_ckernel1);
  cudaMemcpy((void *)&h_ckernel1, (void *)d_ckernel1, sizeof(void *), cudaMemcpyDeviceToHost);
  Iterate(h_ckernel1, 350, 1, d_data);
  cudaDeviceSynchronize();
  cudaFree(d_ckernel1);
  cudaFree(d_data);
  return;
}

int main()
{

  run_test(1);
  run_test(2.0f);

  return 0;
}
