#ifndef __PROTO_ACCEL_HEADER__
#define __PROTO_ACCEL_HEADER__


#ifdef PROTO_ACCEL

#define ACCEL_DECORATION __host__ __device__
#define ACCEL_KERNEL __global__

#else


#define ACCEL_DECORATION /*__host__ __device__ */ 
#define ACCEL_KERNEL /* __global__ */

#endif


#endif
