#include <Proto_gpu.H> //useless

#include <hip/hcc_detail/host_defines.h>
#include <hip/hcc_detail/math_functions.h>
#include <hip/hcc_detail/device_functions.h>

#include<cstdio>
#include<iostream>
#include<string>

template<typename T, typename... Args>
inline void printCPU(T message, Args... args)
{
	printCPU(message);
	printCPU(args...);
}

template<typename T>
inline void printCPU(T message)
{
	std::cout << " CPU: " << message << std::endl;
}

#define FULL_MASK 0xffffffff
__device__
inline int wrapSumReduce( int value) 
{
	for (int offset = warpSize/2; offset > 0; offset /= 2)
    		//__shfl_down(FULL_MASK, value, offset);
    		value+=__shfl_down(value, offset);
	return value;
}
__global__
void initZero(int *ptr)
{
	*ptr = 0;
}

__global__
void reduceSum(int *ptr, int* out, int N)
{
	int sum = int (0);
	for(int i = threadIdx.x + blockIdx.x*blockDim.x ; i < N ; i+= gridDim.x*blockDim.x )
		sum += ptr[i];

	sum=wrapSumReduce(sum);
	if ((threadIdx.x & (warpSize - 1)) == 0)
    		atomicAdd(out, sum);

//	*out += sum;
}	


int main()
{
	// param
	const unsigned int Size = 64;
	assert(Size > 0);

	printCPU ( "Test: compilation __shfl_down_sync  " );
	printCPU ( "Parameter size:" );
	printCPU (  Size );


	int *val_h = new int[Size];
	int *val_d;
	int *val_out;
	const int rnd = 5;
	for( int i = 0 ; i < Size ; i++ ) val_h[i] = 5;
	
	protoMalloc(&val_d,Size*sizeof(int));
	protoMalloc(&val_out,Size*sizeof(int));

	printCPU ( "Value before ");
	printCPU ( rnd );
	protoMemcpy(val_d, val_h, Size*sizeof(int), protoMemcpyHostToDevice);
	protoLaunchKernel(initZero, 1, 1, val_out);
	protoLaunchKernel(reduceSum, 1, 16, val_d, val_out, Size);
	protoMemcpy(val_h, val_out, Size*sizeof(int), protoMemcpyDeviceToHost);

	int endValue = val_h[0];
	
	if(endValue == Size * rnd)
	{
		printCPU ( " Success " );
		printCPU ( "Value after ");
		printCPU ( endValue );
	}
	else
	{
		printCPU ( " Fail " );
		printCPU ( "Value after " );
		printCPU ( endValue );
	}

	return 0;
}
