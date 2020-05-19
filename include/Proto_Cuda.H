#pragma once

#ifdef HIP
#include "hcc.h"
#include "hip/hip_runtime.h"
#else
#include "cuda.h"
#endif


#ifdef HIP
#define protoStream_t hipStream_t
#else 
#define protoStream_t cudaStream_t
#endif


#ifdef HIP

#define HC(X) HIP_CHECK(X)

#define protoMalloc(PTR,NBYTES)    HC(hipMalloc(&PTR,NBYTES)) 
#define protoFree(PTR)             HC(hipFree(PTR))
#define protoHostAlloc(PTR,NBYTES) HC(hipHostMalloc(&PTR,NBYTES))
#define protoHostFree(PTR)         HC(hipHostFree(PTR))

#define protoMemcpy(to,from,size,copyType)             HC(hipMemcpy(to,from,size,copyType))
#define protoMemcpyAsync(to,from,size,copyType,stream) HC(hipMemcpyAsync(to,from,size,copyType, stream))


#define protoStreamCreate(X)      HC(hipStreamCreate(X))
#define protoStreamDestroy(X)     HC(hipStreamDestroy(X))
#define protoStreamSynchronize()  HC(hipDeviceSynchronize())
#define protoSetDevice(X)         HC(hipSetDevice(X))
#define protoStreamSynchronize(X) HC(hipSetDevice(X))

#else

#define protoHostAlloc(PTR,NBYTES) cudaMalloc((void**)&PTR,NBYTES)
#define protoHostFree(PTR)         cudaFree(PTR)
#define protoHostAlloc(PTR,NBYTES) cudaHostAlloc((void**)&PTR,NBYTES) 
#define protoHostFree(PTR)         cudaHostFree(PTR)

#define protoMemcpy(to,from,size,copyType)             cudaMemcpy(to,from,size,copyType)
#define protoMemcpyAsync(to,from,size,copyType,stream) cudaMemcpyAsync(to,from,size,copyType, stream)

#define protoStreamCreate(X)      cudaStreamCreate(X)
#define protoStreamDestroy(X)     cudaStreamDestroy(X)
#define protoStreamSynchronize()  cudaDeviceSynchronize()
#define protoSetDevice(X)         cudaSetDevice(X)
#define protoStreamSynchronize(X) cudaStreamSynchronize

#endif


/*template <class Kernel, class T...>
void protoLaunchKernel( Kernel Ker, int nbBlocs, int nbThreads, T...& args)
{
	protoLaunchKernel(
			Ker,
			nbBlocs,
			nbThreads,
			0,
			0,
			args...
		);
}

template <class Kernel, class T...>
void protoLaunchKernel( Kernel Ker, int nbBlocs, int nbThreads, size_t smem, T..& args)
{
        protoLaunchKernel(
                        Ker,
                        nbBlocs,
                        nbThreads,
                        smem,
                        0,
                        args...
                );
}*/

template <class Kernel, class T...>
void protoLaunchKernel( Kernel Ker, int nbBlocs, int nbThreads, size_t smem = 0, protoStream_t stream = 0, T...& args)
{
#ifdef HIP
	HIP_CHECK(
                        hipLaunchKernelGGL(
                                Ker,
                                dim3( nbBlocs, 0, 0),
                                dim3( nbThreads, 0, 0),
                                smem,
                                stream,
                                args...
                                );
                 );
#else
        Ker<<<nbBlocs, nbThreads,smem,stream>>>(args...);
#endif
}

