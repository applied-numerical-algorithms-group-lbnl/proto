#pragma once
#include <type_traits>

#define cpuMalloc(ptr,nbBytes)     ptr = (decltype(ptr)) malloc(nbBytes)

#define protoMallocCPU(PTR,NBYTES)    cpuMalloc(PTR, NBYTES)
#define protoFreeCPU(PTR)             free(PTR)

#define protoHostAlloc(a,b)        cpuMalloc(PTR, NBYTES)
#define protoFreeHost(PTR)         free(PTR)
// trick

static bool isGpu=false;
#ifndef PROTO_CUDA
#define protoMemcpyDeviceToDevice isGpu
#define protoMemcpyHostToDevice   isGpu
#define protoMemcpyDeviceToHost   isGpu
#endif
// end trick


// COPY


#define protoMemcpyCPU(to,from,size,copyType)             std::memcpy(to,from,size)
#define protoMemcpyAsyncCPU(to,from,size,copyType,stream) std::memcpy(to,from,size)
#define protoMemcpyFromSymbolCPU(a,b,c,d,e)               std::memcpy(a,b,c)
#define protoMemcpyToSymbolCPU(a,b,c,d,e)                 std::memcpy(a,b,c)

/*
#define protoMemcpy(to,from,size,copyType)             std::memcpy(to,from,size)
#define protoMemcpyAsync(to,from,size,copyType,stream) std::memcpy(to,from,size)
#define protoMemcpyFromSymbol(a,b,c,d,e)               std::memcpy(a,b,c)
#define protoMemcpyToSymbol(a,b,c,d,e)                 std::memcpy(a,b,c)
*/


#define protoDeviceSynchronizeCPU()


//#define protoGetCurrentStream 0

/*
class dim3 { 
	public:
	size_t x; 
	size_t y; 
	size_t z; 
	dim3(size_t a_x,size_t a_y, size_t a_z) : x(a_x),y(a_y),z(a_z) {};
        dim3():x(0),y(0),z(0){};
};
*/


// disable commad
