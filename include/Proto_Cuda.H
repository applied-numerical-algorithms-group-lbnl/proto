#pragma once

#ifdef PROTO_HIP
#include "hip/hip_runtime.h"
#else
#include "cuda.h"
#endif


#ifdef PROTO_HIP
#define protoStream_t hipStream_t
#define protoMemcpyDeviceToDevice hipMemcpyDeviceToDevice
#define protoMemcpyHostToDevice hipMemcpyHostToDevice
#define protoMemcpyDeviceToHost hipMemcpyDeviceToHost
#define protoError hipError_t
#define protoSuccess hipSuccess
#define protoDeviceProp hipDeviceProp_t

#define protoGetLastError()    hipGetLastError()	  
#define protoGetErrorString(X) hipGetErrorString(X)

#ifndef NDEBUG
	#define HIP_CHECK(in)                                    \
    	do                                                          \
    	{                                                           \
        	hipError_t error = in;                           \
        	if(error == hipSuccess) std::cout << hipGetErrorString(error); \
		exit(0); \
    	} while(0)
#else
        #define HIP_CHECK(condition) condition
#endif

#else 
#define protoStream_t cudaStream_t
#define protoMemcpyDeviceToDevice cudaMemcpyDeviceToDevice
#define protoMemcpyHostToDevice cudaMemcpyHostToDevice
#define protoMemcpyDeviceToHost cudaMemcpyDeviceToHost
#define protoError cudaError
#define protoSuccess cudaSuccess
#define protoDeviceProp cudaDeviceProp

#endif

#ifdef PROTO_HIP

#define HC(X) HIP_CHECK(X)

#define protoMalloc(PTR,NBYTES)    HC(hipMalloc(PTR,NBYTES)) 
#define protoFree(PTR)             HC(hipFree(PTR))
#define protoHostAlloc(PTR,NBYTES) HC(hipHostMalloc(PTR,NBYTES))
#define protoFreeHost(PTR)         HC(hipHostFree(PTR))
#define protoHostFree(PTR)         protoFreeHost(PTR)

#define protoMemcpy(to,from,size,copyType)             HC(hipMemcpy(to,from,size,copyType))
#define protoMemcpyAsync(to,from,size,copyType,stream) HC(hipMemcpyAsync(to,from,size,copyType, stream))
#define protoMemcpyFromSymbol(a,b,c,d,e)               hipMemcpyFromSymbol(a,b,c,d,e)


#define protoStreamCreate(X)          HC(hipStreamCreate(X))
#define protoStreamDestroy(X)         HC(hipStreamDestroy(X))
#define protoDeviceSynchronize()      HC(hipDeviceSynchronize())
#define protoStreamSynchronize(X)     HC(hipStreamSynchronize(X))
#define protoSetDevice(X)             HC(hipSetDevice(X))
#define protoGetDeviceProperties(X,Y) HC(hipGetDeviceProperties(X,Y))
#else

#define protoMalloc(a,b)           cudaMalloc(a,b)
#define protoFree(a)               cudaFree(a)
#define protoHostAlloc(a,b)        cudaMallocHost(a,b) 
#define protoFreeHost(PTR)         cudaFreeHost(PTR)

#define protoMemcpy(to,from,size,copyType)             cudaMemcpy(to,from,size,copyType)
#define protoMemcpyAsync(to,from,size,copyType,stream) cudaMemcpyAsync(to,from,size,copyType, stream)
#define protoMemcpyFromSymbol(a,b,c,d,e)               cudaMemcpyFromSymbol(a,b,c,d,e)

#define protoGetLastError()    cudaGetLastError()	  
#define protoGetErrorString(X) cudaGetErrorString(X)

#define protoStreamCreate(X)          cudaStreamCreate(X)
#define protoStreamDestroy(X)         cudaStreamDestroy(X)
#define protoDeviceSynchronize()      cudaDeviceSynchronize()
#define protoStreamSynchronize(X)     cudaStreamSynchronize(X)
#define protoSetDevice(X)             cudaSetDevice(X)
#define protoGetDeviceProperties(X,Y) cudaGetDeviceProperties(X,Y)

#endif


template <class Kernel, typename... T>
void protoLaunchKernel( Kernel Ker, int nbBlocks, int nbThreads, T... args)
{
#ifdef PROTO_HIP
	hipLaunchKernelGGL(
		Ker,
		dim3( nbBlocks, 0, 0),
		dim3( nbThreads, 0, 0),
		0,
		0,
		args...
	);
#else
        Ker<<<nbBlocks, nbThreads,0,0>>>(args...);
#endif
}

template <class Kernel, typename... T>
void protoLaunchKernelMem( Kernel Ker, int nbBlocks, int nbThreads, size_t smem, T... args)
{
#ifdef PROTO_HIP
	hipLaunchKernelGGL(
		Ker,
		dim3( nbBlocks, 0, 0),
		dim3( nbThreads, 0, 0),
		smem,
		0,
		args...
        );
#else
        Ker<<<nbBlocks, nbThreads,smem,0>>>(args...);
#endif
}

template <class Kernel, typename... T>
void protoLaunchKernelMemAsync( Kernel Ker, size_t nbBlocks, size_t nbThreads, size_t smem, protoStream_t stream, T... args)
{
#ifdef PROTO_HIP
	hipLaunchKernelGGL(
		Ker, 
		dim3( nbBlocks, 0, 0), 
		dim3( nbThreads, 0, 0), 
		smem, 
		stream, 
		args...
	);
#else
        Ker<<<nbBlocks, nbThreads,smem,stream>>>(args...);
#endif
}
