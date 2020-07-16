/// In general, all calls to cuda have been redefined as macros in order to easily make a switch with hip at compile time. 
/// These macros are defined into \texttt{Proto\_gpu.H}. 
/// To use the HIP, you need to include the following flags: \texttt{-DPROTO\_CUDA} and \texttt{-DPROTO\_HIP}.
/// Most of cudaNAME functions are renamed protoNAME such as:
///
/// \begin{lstlisting}[language=C++,caption={Macro Define}]
/// #ifdef PROTO_HIP
/// #define protoMalloc(...) hipMalloc(...) // HIP
/// #else   protoMalloc(...) cudaMalloc(...) // CUDA
/// #endif
/// \end{lstlisting}
///
/// In the following sections, you will find the renamed cuda functions and data types. 
/// Functions that aren't real cuda Function such as cudaApply or structures such as cudaUglyStruct keep their names. 
/// More than 95$\%$ of the changes have been made in the Proto.


#pragma once

#ifdef PROTO_CUDA


#ifdef PROTO_HIP
#include "hip/hip_runtime.h"
#include "hip/hip_runtime_api.h"
#include "iostream" // Use for CHECK(X)
#else
#include "cuda.h"
#endif


/// Data Types / classes
#ifdef PROTO_HIP // HIP
	#define protoStream_t             hipStream_t
	#define protoMemcpyDeviceToDevice hipMemcpyDeviceToDevice
	#define protoMemcpyHostToDevice   hipMemcpyHostToDevice
	#define protoMemcpyDeviceToHost   hipMemcpyDeviceToHost
	#define protoError                hipError_t //legacy
	#define protoError_t              hipError_t
	#define protoSuccess              hipSuccess
	#define protoDeviceProp           hipDeviceProp_t
	#define protoPointerAttributes    hipPointerAttribute_t
	#define protoPitchedPtr           hipPitchedPtr
	#define protoArray                hipArray
	#define protoExtent               hipExtent
	#define protoChannelFormatDesc    hipChannelFormatDesc
	#define protoReadModeElementType  hipReadModeElementType


	#define protoGetLastError()    hipGetLastError()	  
	#define protoGetErrorString(X) hipGetErrorString(X)

	#ifndef NDEBUG
		#define HIP_CHECK(in)                                    \
		do                                                          \
		{                                                           \
			hipError_t error = in;                           \
			if(error != hipSuccess) \
			{ \
				std::cout << hipGetErrorString(error); \
				exit(0); \
			}\
		} while(0)
	#else
		#define HIP_CHECK(condition) condition
	#endif

#else // CUDA
	#define protoStream_t             cudaStream_t
	#define protoMemcpyDeviceToDevice cudaMemcpyDeviceToDevice
	#define protoMemcpyHostToDevice   cudaMemcpyHostToDevice
	#define protoMemcpyDeviceToHost   cudaMemcpyDeviceToHost
	#define protoError                cudaError // legacy
	#define protoError_t              cudaError
	#define protoSuccess              cudaSuccess
	#define protoDeviceProp           cudaDeviceProp
	#define protoPointerAttributes    cudaPointerAttributes
	#define protoPitchedPtr           cudaPitchedPtr
	#define protoArray                cudaArray
	#define protoExtent               cudaExtent
	#define protoChannelFormatDesc    cudaChannelFormatDesc
	#define protoReadModeElementType  cudaReadModeElementType
#endif


/// Functions
#ifdef PROTO_HIP // HIP

	#define HC(X) HIP_CHECK(X)

	// MEMORY
	#define protoMalloc(PTR,NBYTES)    HC(hipMalloc(PTR,NBYTES)) 
	#define protoFree(PTR)             HC(hipFree(PTR))
	#define protoHostAlloc(PTR,NBYTES) HC(hipHostMalloc(PTR,NBYTES))
	#define protoFreeHost(PTR)         HC(hipHostFree(PTR))
	#define protoHostFree(PTR)         hipFreeHost(PTR)
	#define protoMallocManaged(a,b)    HC(hipMallocManaged(a,b))
	#define protoMemset(a,b,c)         HC(hipMemset(a,b,c))

	// COPY
	#define protoMemcpy(to,from,size,copyType)             HC(hipMemcpy(to,from,size,copyType))
	#define protoMemcpyAsync(to,from,size,copyType,stream) HC(hipMemcpyAsync(to,from,size,copyType, stream))
	#define protoMemcpyFromSymbol(a,b,c,d,e)               hipMemcpyFromSymbol(a,b,c,d,e)
	#define protoMemcpyToSymbol(a,b,c,d,e)                 hipMemcpyToSymbol(a,b,c,d,e)

	// STREAM
	#define protoStreamCreate(X)           HC(hipStreamCreate(X))
	#define protoStreamDestroy(X)          HC(hipStreamDestroy(X))
	#define protoDeviceSynchronize()       hipDeviceSynchronize()
	#define protoStreamSynchronize(X)      HC(hipStreamSynchronize(X))

	// DEVICE
	#define protoSetDevice(X)              HC(hipSetDevice(X))
	#define protoGetDeviceProperties(X,Y)  HC(hipGetDeviceProperties(X,Y))
	#define protoDeviceReset()             HC(hipDeviceReset())
	#define protoPointerGetAttributes(X,Y) HC(hipPointerGetAttributes(X,Y))

	// OTHER
	#define protoBindTexture(a,b,c,d,e)    HC(hipBindTexture(a,b,c,d,e))
	#define protoMalloc3D(a,b) HC(hipMalloc3D(a,b))
	#define make_protoExtent   hip_cudaExtent

#else //CUDA
	// MEMORY
	#define protoMalloc(a,b)           cudaMalloc(a,b)
	#define protoFree(a)               cudaFree(a)
	#define protoHostAlloc(a,b)        cudaMallocHost(a,b) 
	#define protoFreeHost(PTR)         cudaFreeHost(PTR)
	#define protoMallocManaged(a,b)    cudaMallocManaged(a,b)
	#define protoMemset(a,b,c)         cudaMemset(a,b,c)

	// COPY
	#define protoMemcpy(to,from,size,copyType)             cudaMemcpy(to,from,size,copyType)
	#define protoMemcpyAsync(to,from,size,copyType,stream) cudaMemcpyAsync(to,from,size,copyType, stream)
	#define protoMemcpyFromSymbol(a,b,c,d,e)               cudaMemcpyFromSymbol(a,b,c,d,e)
	#define protoMemcpyToSymbol(a,b,c,d,e)                 cudaMemcpyToSymbol(a,b,c,d,e)

	// ERROR
	#define protoGetLastError()    cudaGetLastError()	  
	#define protoGetErrorString(X) cudaGetErrorString(X)

	// STREAM
	#define protoStreamCreate(X)           cudaStreamCreate(X)
	#define protoStreamDestroy(X)          cudaStreamDestroy(X)
	#define protoDeviceSynchronize()       cudaDeviceSynchronize()
	#define protoStreamSynchronize(X)      cudaStreamSynchronize(X)

	// DEVICE
	#define protoSetDevice(X)              cudaSetDevice(X)
	#define protoGetDeviceProperties(X,Y)  cudaGetDeviceProperties(X,Y)
	#define protoDeviceReset()             cudaDeviceReset()
	#define protoPointerGetAttributes(X,Y) cudaPointerGetAttributes(X,Y)

	// OTHER
	#define protoBindTexture(a,b,c,d,e) cudaBindTexture(a,b,c,d,e)
	#define protoMalloc3D(a,b) cudaMalloc3D(a,b)
	#define make_protoExtent   make_cudaExtent
#endif

/*
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
}*/

#ifdef PROTO_HIP
#define protoLaunchKernel(Ker, nbBlocks, nbThreads, args...) hipLaunchKernelGGL( \
                Ker, \
                nbBlocks, \
                nbThreads, \
                0, \
                0, \
                args)
#else
#define protoLaunchKernel(Ker, nbBlocks, nbThreads, args...) Ker<<<nbBlocks,nbThreads>>>(args)
#endif

#ifdef PROTO_HIP
#define protoLaunchKernelMem(Ker, nbBlocks, nbThreads, smem, stream, args...) hipLaunchKernelGGL( \
                Ker, \
                dim3( nbBlocks, 0, 0), \
                dim3( nbThreads, 0, 0), \
                smem, \
                0, \
                args)
#else
#define protoLaunchKernelMem(Ker, nbBlocks, nbThreads, smem, args...) Ker<<<nbBlocks, nbThreads,smem>>(args)
#endif



#ifdef PROTO_HIP
#define protoLaunchKernelMemAsync(Ker, nbBlocks, nbThreads, smem, stream, args...) hipLaunchKernelGGL( \
                Ker, \
                nbBlocks, \
                nbThreads, \
                smem, \
                stream, \
                args)
#else
#define protoLaunchKernelMemAsync(Ker, nbBlocks, nbThreads, smem, stream, args...) Ker<<<nbBlocks, nbThreads,smem,stream>>>(args)
#endif

#endif
