#include<Proto_gpu.H>
#include<iostream>

#define PROTO_KERNEL_START inline __device__

#define PROTO_KERNEL_END(local_name, app_name)                  \
 typedef struct { \
    const char* myname = #app_name; \
    template <typename... T> \
    inline __device__ void operator()(T... args){ local_name(args...);}; \
  } struct_##local_name; \
  static struct_##local_name app_name;

PROTO_KERNEL_START
void myAddFunctionF(double* array, size_t idx, double value)
        {
                array[idx]+=value;
        }
PROTO_KERNEL_END(myAddFunctionF, myAddFunction)

template<typename Struct, typename T, typename... Args>
__global__
void meshStructLauncher(Struct st, T* arr, size_t end, Args... args)
{
        int idx = threadIdx.x + blockIdx.x * blockDim.x ;

        if(idx<end)
                st(arr,idx,args...);
}

template< typename Struct, typename T, typename... Args>
void structLauncher(Struct &st, T* in, int end, Args&&... args)
{

        T* __restrict__ data = in;

        const int threads = 256;
        const int blocks = end / threads + 1;

        protoLaunchKernel(meshStructLauncher, blocks, threads, st, data, end, std::forward<Args>(args)...);
}


//#define comparaison
#ifdef comparaison
__global__
void myAddFunctionNoForall(double* array, size_t size, double value)
{
	const unsigned int tid = threadIdx.x + blockDim.x*blockIdx.x;
	if(tid<size)
		array[tid]+=value;
}
#endif

int main()
{
	int size = 10000;
	// allocations
	double *hostPtr = new double[size];
	double *devicePtr;
	protoMalloc(&devicePtr, size * sizeof(double));

	// init on host
	for(int i = 0 ; i<size ; i++)
		hostPtr[i] = 3;

	protoMemcpy(devicePtr, hostPtr, size * sizeof(double) , protoMemcpyHostToDevice);

#ifdef comparaison // cuda
	std::cout << " We are using the basic add function without the forall design " <<std::endl;
	myAddFunctionNoForall<<<(size+256-1)/256, 256>>>(devicePtr, size, 663);
#else
	structLauncher( myAddFunction, devicePtr, size, 663);
#endif

	protoMemcpy(hostPtr, devicePtr, size * sizeof(double) , protoMemcpyDeviceToHost);


	for(int i = 0 ; i<size ; i++)
		if(hostPtr[i] != 666 )
		       {
		       	       std::cout << " Failled: [" << i << "] " << hostPtr[i] << std::endl;
				break;
		       }

	std::cout << " fin " << std::endl;
	return 0;
}

