#ifndef __PROTO_MACROS_H_
#define __PROTO_MACROS_H_

#include "Proto_MemType.H"
#ifdef PROTO_CUDA
#include "Proto_gpu.H"
#endif
#include "Proto_cpu.H"

#ifndef PROTO_CUDA
//Some of the device-agnostic code needs definitions of dim3, protoStream_t, protoGetCurrentStream
class dim3 { 
	public:
	size_t x; 
	size_t y; 
	size_t z; 
	dim3(size_t a_x,size_t a_y, size_t a_z) : x(a_x),y(a_y),z(a_z) {};
};
#define protoStream_t int
#define protoGetCurrentStream 0

#define protoMalloc(memtype,PTR,NBYTES) protoMallocCPU(PTR,NBYTES)
#define protoFree(memtype,PTR) protoFreeCPU(PTR)

#define protoMemcpy(memtype,to,from,size,copyType) protoMemcpyCPU(to,from,size,copyType)
#define protoMemcpyAsync(memtype,to,from,size,copyType,stream) protoMemcpyAsyncCPU(to,from,size,copyType,stream)
#define protoMemcpyFromSymbol(memtype,a,b,c,d,e) protoMemcpyFromSymbolCPU(a,b,c,d,e)
#define protoMemcpyToSymbol(memtype,a,b,c,d,e) protoMemcpyToSymbolCPU(a,b,c,d,e)

#define protoDeviceSynchronize(memtype) protoDeviceSynchronizeCPU()

template<Proto::MemType Type, typename Ker, typename... Args>
inline void protoLaunchKernelT(int nBlocks, int nbThreads, const Args&... args)
{
  protoLaunchKernelCPU(Ker::cpu, nBlock, nbThreads, args...);
}

template<Proto::MemType Type, typename Ker, typename... Args>
inline void protoLaunchKernelT(dim3 nBlocks, int nbThreads, const Args&... args)
{
  protoLaunchKernelCPU(Ker::cpu, nBlock, nbThreads, args...);
}

//#define protoLaunchKernel(memtype,Ker, nbBlocks, nbThreads, args...)  \
//    protoLaunchKernelCPU(Ker, nbBlocks, nbThreads, args)
#define protoLaunchKernelMem(memtype,Ker, nbBlocks, nbThreads, smem, args...) \
    protoLaunchKernelMemCPU(Ker, nbBlocks, nbThreads, smem, args)

template<Proto::MemType Type, typename Ker, typename... Args>
inline void protoLaunchKernelMemAsyncT(int nBlocks, int nbThreads, int smem, protoStream_t stream, const Args&... args)
{
  protoLaunchKernelMemAsyncCPU(Ker::cpu, nbBlocks, nbThreads, smem, stream, args...);
}

template<Proto::MemType Type, typename Ker, typename... Args>
inline void protoLaunchKernelMemAsyncT(dim3 nBlocks, int nbThreads, int smem, protoStream_t stream, const Args&... args)
{
  protoLaunchKernelMemAsyncCPU(Ker::cpu, nbBlocks, nbThreads, smem, stream, args...);
}

#define protoLaunchKernelMemAsync(memtype,Ker, nbBlocks, nbThreads, smem, stream, args...) \
    protoLaunchKernelMemAsyncCPU(Ker, nbBlocks, nbThreads, smem, stream, args)

#else //If PROTO_CUDA is defined

#define protoMalloc(memtype,PTR,NBYTES) \
    if(memtype==Proto::MemType::DEVICE) \
      protoMallocGPU(PTR,NBYTES);	\
    else \
      protoMallocCPU(PTR,NBYTES);

#define protoFree(memtype,PTR) \
    if(memtype==Proto::MemType::DEVICE) \
      protoFreeGPU(PTR);		\
    else \
      protoFreeCPU(PTR);

#define protoMemcpy(memtype,to,from,size,copyType) \
    if(memtype==Proto::MemType::DEVICE) \
      protoMemcpyGPU(to,from,size,copyType);	\
    else \
      protoMemcpyCPU(to,from,size,copyType);
#define protoMemcpyAsync(memtype,to,from,size,copyType,stream) \
    if(memtype==Proto::MemType::DEVICE) \
      protoMemcpyAsyncGPU(to,from,size,copyType,stream);	\
    else \
      protoMemcpyAsyncCPU(to,from,size,copyType,stream);
#define protoMemcpyFromSymbol(memtype,a,b,c,d,e) \
    if(memtype==Proto::MemType::DEVICE) \
      protoMemcpyFromSymbolGPU(a,b,c,d,e);	\
    else \
      protoMemcpyFromSymbolCPU(a,b,c,d,e);
#define protoMemcpyToSymbol(memtype,a,b,c,d,e) \
    if(memtype==Proto::MemType::DEVICE) \
      protoMemcpyToSymbolGPU(a,b,c,d,e);	\
    else \
      protoMemcpyToSymbolCPU(a,b,c,d,e);

#define protoDeviceSynchronize(memtype) \
    if(memtype==Proto::MemType::DEVICE) \
      protoDeviceSynchronizeGPU();	\
    else \
      protoDeviceSynchronizeCPU();

#define protoLaunchKernel(memtype,Ker, nbBlocks, nbThreads, args...) \
    if(memtype==Proto::MemType::DEVICE) \
        protoLaunchKernelGPU(Ker##GPU, nbBlocks, nbThreads, args) \
    else \
        protoLaunchKernelCPU(Ker##CPU, nbBlocks, nbThreads, args)
#define protoLaunchKernelMem(memtype,Ker, nbBlocks, nbThreads, smem, args...) \
    if(memtype==Proto::MemType::DEVICE) \
      protoLaunchKernelMemGPU(Ker##GPU, nbBlocks, nbThreads, smem, args);	\
    else \
      protoLaunchKernelMemCPU(Ker##CPU, nbBlocks, nbThreads, smem, args);
#define protoLaunchKernelMemAsync(memtype,Ker, nbBlocks, nbThreads, smem, stream, args...) \
    if(memtype==Proto::MemType::DEVICE) \
      protoLaunchKernelMemAsyncGPU(Ker##GPU, nbBlocks, nbThreads, smem, stream, args); \
    else \
      protoLaunchKernelMemAsyncCPU(Ker##CPU, nbBlocks, nbThreads, smem, stream, args);

#endif

/*
#define GPU(name) __global__ \
void gpu_##name

#define CPU(name) void cpu_##name

#define FUNCTOR(name) template<> struct base##name{     \
    };\
\
    template<> struct base##name<true>\
    {\
    template<typename... T>\
    void operator()(unsigned int nbBlocks, unsigned int nbThreads, T... args)\
    {\
    gpu_##name<<<nbBlocks,nbThreads>>>(args...);\
    }\
    };\
    template<> struct base##name<false>\
    {\
    template<typename... T>\
    void operator()(unsigned int nbBlocks, unsigned int nbThreads, T... args)\
    {\
    cpu_##name(args...);\
    }\
    };\
    base##name<true> gpu##name;\
    base##name<false> cpu##name;
*/

#endif
