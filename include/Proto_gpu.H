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

#define protoGetCurrentStream DisjointBoxLayout::getCurrentStream()

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
	#define protoEvent_t		  hipEvent_t 
	#define protoGetLastError()       hipGetLastError()	  
	#define protoPeekAtLastError()    hipPeekAtLastError()
	#define protoGetErrorString(X)    hipGetErrorString(X)
	#define protoThreadSynchronize()  hipDeviceSynchronize()


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
	#define protoEvent_t		  cudaEvent_t 
	#define protoGetLastError()       cudaGetLastError()	  
	#define protoPeekAtLastError()    cudaPeekAtLastError()
	#define protoGetErrorString(X)    cudaGetErrorString(X)
	#define protoThreadSynchronize()  cudaThreadSynchronize()

#endif

#ifndef NDEBUG

	#ifndef superDebug
		#define GPU_CHECK(in)                                    \
		do                                                          \
		{                                                           \
			protoError_t error = in;                           \
			if(error != protoSuccess) \
			{ \
				std::cout << protoGetErrorString(error); \
				exit(0); \
			}\
		} while(0)
	#else
		#define GPU_CHECK(in)                                    \
                do                                                          \
                {                                                           \
			std::cout << "Try " << #in << " file: "<< __FILE__ << " line: " << __LINE__  << std::endl; \
                        protoError_t error = in;                           \
			protoDeviceSynchronizeGPU();\
                        if(error != protoSuccess) \
                        { \
                                std::cout << protoGetErrorString(error); \
                                exit(0); \
                        }\
			else std::cout << "Success " << #in << std::endl; \
                } while(0)
	#endif

#else
	#define GPU_CHECK(condition) condition
#endif

/// Functions
#ifdef PROTO_HIP // HIP

	#define HC(X) GPU_CHECK(X)

	// MEMORY
	#define protoMallocGPU(PTR,NBYTES)    HC(hipMalloc(&PTR,NBYTES)) 
	#define protoFreeGPU(PTR)             HC(hipFree(PTR))
	#define protoHostAlloc(PTR,NBYTES) HC(hipHostMalloc(&PTR,NBYTES))
	#define protoFreeHost(PTR)         HC(hipHostFree(PTR))
	#define protoHostFree(PTR)         HC(hipFreeHost(PTR))
	#define protoMallocManaged(a,b)    HC(hipMallocManaged(&a,b))
	#define protoMemset(a,b,c)         HC(hipMemset(a,b,c))

	// COPY
	#define protoMemcpyGPU(to,from,size,copyType)             HC(hipMemcpy(to,from,size,copyType))
	#define protoMemcpyAsyncGPU(to,from,size,copyType,stream) HC(hipMemcpyAsync(to,from,size,copyType, stream))
	#define protoMemcpyFromSymbolGPU(a,b,c,d,e)               hipMemcpyFromSymbol(a,b,c,d,e) // not used anymore
	#define protoMemcpyToSymbolGPU(a,b,c,d,e)                 hipMemcpyToSymbol(a,b,c,d,e)

	// STREAM
#define protoDeviceSynchronizeGPU()       hipDeviceSynchronize()
#define protoStreamCreate(X)           HC(hipStreamCreate(X))
	#define protoStreamDestroy(X)          HC(hipStreamDestroy(X))
	#define protoStreamSynchronize(X)      HC(hipStreamSynchronize(X))

	// DEVICE
	#define protoSetDevice(X)              HC(hipSetDevice(X))
	#define protoGetDeviceProperties(X,Y)  HC(hipGetDeviceProperties(X,Y))
	#define protoDeviceReset()             HC(hipDeviceReset())
	#define protoPointerGetAttributes(X,Y) HC(hipPointerGetAttributes(X,Y))
	#define protoGetDeviceCount(X)         HC(hipGetDeviceCount(X))
	#define protoGetDevice(X)              HC(hipGetDevice(X))

	// EVENT
	#define protoEventCreate(X) 		HC(hipEventCreate(X))
	#define protoEventRecord(X) 		HC(hipEventRecord(X))
	#define protoEventSynchronize(X) 	hipEventSynchronize(X)
	#define protoEventElapsedTime(a,b,c) 	HC(hipEventElapsedTime(a,b,c))

	// OTHER
	#define protoBindTexture(a,b,c,d,e)    HC(hipBindTexture(a,b,c,d,e))
	#define protoMalloc3D(a,b) HC(hipMalloc3D(&a,b))
	#define make_protoExtent   hip_cudaExtent

#else //CUDA

	#define CC(X) GPU_CHECK(X) // CudaCheck

	// MEMORY
	#define protoMallocGPU(a,b)           CC(cudaMalloc(&a,b))
	#define protoFreeGPU(a)               CC(cudaFree(a))
	#define protoHostAlloc(a,b)        CC(cudaMallocHost(&a,b))
	#define protoFreeHost(PTR)         CC(cudaFreeHost(PTR))
	#define protoMallocManaged(a,b)    CC(cudaMallocManaged(&a,b))
	#define protoMemset(a,b,c)         CC(cudaMemset(a,b,c))

	// COPY
	#define protoMemcpyGPU(to,from,size,copyType)             CC(cudaMemcpy(to,from,size,copyType))
	#define protoMemcpyAsyncGPU(to,from,size,copyType,stream) CC(cudaMemcpyAsync(to,from,size,copyType, stream))
	#define protoMemcpyFromSymbolGPU(a,b,c,d,e)               cudaMemcpyFromSymbol(a,b,c,d,e)
	#define protoMemcpyToSymbolGPU(a,b,c,d,e)                 cudaMemcpyToSymbol(a,b,c,d,e)

	// STREAM
	#define protoDeviceSynchronizeGPU()       cudaDeviceSynchronize()
	#define protoStreamCreate(X)           CC(cudaStreamCreate(X))
	#define protoStreamDestroy(X)          CC(cudaStreamDestroy(X))
	#define protoStreamSynchronize(X)      CC(cudaStreamSynchronize(X))

	// DEVICE
	#define protoSetDevice(X)              CC(cudaSetDevice(X))
	#define protoGetDeviceProperties(X,Y)  CC(cudaGetDeviceProperties(X,Y))
	#define protoDeviceReset()             CC(cudaDeviceReset())
	#define protoPointerGetAttributes(X,Y) CC(cudaPointerGetAttributes(X,Y))
	#define protoGetDevice(X)              CC(cudaGetDevice(X))
	#define protoGetDeviceCount(X)         CC(cudaGetDeviceCount(X))

	// EVENT
	#define protoEventCreate(X) 		CC(cudaEventCreate(X))
	#define protoEventRecord(X) 		CC(cudaEventRecord(X))
	#define protoEventSynchronize(X) 	cudaEventSynchronize(X)
	#define protoEventElapsedTime(a,b,c) 	CC(cudaEventElapsedTime(a,b,c))

	// OTHER
	#define protoBindTexture(a,b,c,d,e) CC(cudaBindTexture(a,b,c,d,e))
	#define protoMalloc3D(a,b) CC(cudaMalloc3D(&a,b))
	#define make_protoExtent   make_cudaExtent
#endif

#include <iostream>

#ifdef superDebug
#define PRINT_KER(X) std::cout << "Kernel: "<< #X  << " file " << __FILE__ << " line " <<__LINE__<< std::endl; \
		X \
		protoDeviceSynchronizeGPU();\
		{protoError_t error = protoPeekAtLastError();                    \
		protoDeviceSynchronizeGPU();\
                if(error != protoSuccess) \
                { \
                  std::cout << protoGetErrorString(error); \
                  exit(0); \
                }\
		else std::cout << "Success Kernel: "<< #X << std::endl;}
#define PRINT_KER_CUDA(X) std::cout << "Kernel: "<< #X  << " file " << __FILE__ << " line " <<__LINE__<< std::endl; \
		X; \
		protoDeviceSynchronizeGPU();\
		{protoError_t error = protoPeekAtLastError();                    \
		protoDeviceSynchronizeGPU();\
                if(error != protoSuccess) \
                { \
                  std::cout << protoGetErrorString(error); \
                  exit(0); \
                }\
		else std::cout << "Success Kernel: "<< #X << std::endl;}
#else
#define PRINT_KER(X) X
#define PRINT_KER_CUDA(X) X
#endif

#ifdef PROTO_HIP
#define protoLaunchKernelGPU(Ker, nbBlocks, nbThreads, args...) PRINT_KER(hipLaunchKernelGGL(Ker, nbBlocks, nbThreads, 0, 0, args);)
#else
#define protoLaunchKernelGPU(Ker, nbBlocks, nbThreads, args...) PRINT_KER_CUDA( ( Ker<<<nbBlocks,nbThreads>>>(args)) )
#endif

#ifdef PROTO_HIP
#define protoLaunchKernelMemGPU(Ker, nbBlocks, nbThreads, smem, args...) PRINT_KER(hipLaunchKernelGGL( Ker, nbBlocks, nbThreads, smem, 0, args);)
#else
#define protoLaunchKernelMemGPU(Ker, nbBlocks, nbThreads, smem, args...) PRINT_KER_CUDA(( Ker<<<nbBlocks, nbThreads,smem>>(args)))
#endif



#ifdef PROTO_HIP
#define protoLaunchKernelMemAsyncGPU(Ker, nbBlocks, nbThreads, smem, stream, args...) PRINT_KER( hipLaunchKernelGGL( Ker, nbBlocks, nbThreads, smem, stream, args);) 
#else
#define protoLaunchKernelMemAsyncGPU(Ker, nbBlocks, nbThreads, smem, stream, args...) PRINT_KER_CUDA((Ker<<<nbBlocks, nbThreads,smem,stream>>>(args)))
#endif

#ifdef PROTO_CUDA
template<typename T>
inline bool isDeviceMemory(T* ptr)
{
  protoPointerAttributes att;
  protoPointerGetAttributes(&att, ptr);
#ifdef PROTO_HIP
  if(att.memoryType == hipMemoryTypeDevice) return true; // = 2-> device allocation
#else
  if(att.type == 2) return true; // = 2-> device allocation
#endif
  return false;
}

#endif

//// tuning
// 

inline void v100tuning(int nbElems, int & nbBlocks, int &blockSize)
{
  // determine the best block size and block dim
  blockSize = 256;
  nbBlocks = ( nbElems + blockSize - 1)/ blockSize;
  if(nbBlocks < 80)
  {
      // On V100 we want at least 80 blocks;
      nbBlocks = 80;
      // figure out what is the blockSize for 80 blocks
      // 
      blockSize = ( nbElems + nbBlocks - 1) / nbBlocks;
      // as we use __syncthreads(), we want that stride modulo 32 is equal to 0
      int coeff = blockSize / 32;
      if(coeff > 0)
      {
        blockSize = coeff * 32;
        // recompute the number of blocks > 80
        nbBlocks = ( nbElems + blockSize - 1) / nbBlocks;
      }
  }
}
#endif
