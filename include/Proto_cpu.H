#pragma once
#include <type_traits>

#define cpuMalloc(ptr,nbBytes)     ptr = (decltype(ptr)) malloc(nbBytes)

#define protoMallocCPU(PTR,NBYTES)    cpuMalloc(PTR, NBYTES)
#define protoFreeCPU(PTR)             free(PTR)

// trick

/*static bool isGpu false;

#define protoMemcpyDeviceToDevice isGpu
#define protoMemcpyHostToDevice   isGpu
#define protoMemcpyDeviceToHost   isGpu
*/
// end trick


// COPY
#define protoMemcpyCPU(to,from,size,copyType)             std::memcpy(to,from,size)
#define protoMemcpyAsyncCPU(to,from,size,copyType,stream) std::memcpy(to,from,size)
#define protoMemcpyFromSymbolCPU(a,b,c,d,e)               std::memcpy(a,b,c)
#define protoMemcpyToSymbolCPU(a,b,c,d,e)                 std::memcpy(a,b,c)

#define protoDeviceSynchronizeCPU()

#define protoLaunchKernelCPU(Ker, nbBlocks, nbThreads, args...) Ker(args)
#define protoLaunchKernelMemCPU(Ker, nbBlocks, nbThreads, smem, args...) Ker(args)
#define protoLaunchKernelMemAsyncCPU(Ker, nbBlocks, nbThreads, smem, stream, args...) Ker(args)

// disable commad
