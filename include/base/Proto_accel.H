#ifndef __PROTO_ACCEL_HEADER__
#define __PROTO_ACCEL_HEADER__

#if defined PROTO_ACCEL && !defined PROTO_SYCL

#define ACCEL_DECORATION __host__ __device__
#define ACCEL_KERNEL __global__

#else

#define ACCEL_DECORATION /*__host__ __device__ */ 
#define ACCEL_KERNEL /* __global__ */

#endif

#ifdef PROTO_SYCL
#include <sycl/sycl.hpp>
#endif

#endif
