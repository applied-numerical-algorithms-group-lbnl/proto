#ifndef __PROTO_REDUCTION_H__
#define __PROTO_REDUCTION_H__

#include "Proto_gpu.H"

#include <limits>
#include <typeinfo>

namespace Proto {

enum Operation { Max, Min, Abs };
enum Atomic { Warp, Block, None };

constexpr int line = 128; // bytes in a cache line
template<typename T>
CUDA_DECORATION
static constexpr T max(T a, T b) { return a > b ? a : b; }
template<typename T>
CUDA_DECORATION
static constexpr T min(T a, T b) { return a < b ? a : b; }

template<typename T, Operation op>
CUDA_DECORATION // every block's sum is written to out[]
T compare(T last, T next) {

    T res;
    switch(op) {
        case Max:
            res = max<T>(last,next);
	    break;
        case Min:
            res = min<T>(last,next);
	    break;
        case Abs:
            res = (next > 0 ? max<T>(last,next) : max<T>(last,-next));
	    break;
    }
   return res;
}

template<typename T, Operation op>
CUDA_DECORATION
T reset() {
    T reset;
    switch (op) {
        case Max:
            reset = std::numeric_limits<T>::min();
            break;
        case Min:
            reset = std::numeric_limits<T>::max();
            break;
        case Abs:
            reset = static_cast<T>(0);
            break;
    }
    return reset;
}

#ifdef PROTO_CUDA

template<typename T, Operation op>
__device__ // reduction within a warp, only lane 0 gets right result
T warpOp(T val, int idx, int size) {
    unsigned mask = 0xffffffff; // FULL_MASK
    if ((size/warpSize)*warpSize <= idx) // is thread in the top warp?
        mask = (1 << size%warpSize) - 1; // set bits = # of active threads in top warp
    for (unsigned int delta = warpSize/2; delta > 0; delta /= 2)
#ifdef PROTO_HIP
        val = compare<T,op>(val, (T)__shfl_down(val,delta)); // __shfl_down depreciated since CUDA 9.0 but __shfl_down_sync isn't available with hipcc (17/07/2020)
#else
        val = compare<T,op>(val, (T)__shfl_down_sync(mask,val,delta));
#endif
    return val;
}

template<typename T, Operation op>
__device__ // reduce within a block by storing partial warp reductions in shmem
T blockOp(T val, int idx, int size) {
    extern __shared__ __align__(sizeof(T)) unsigned char shdata[];
    T *shmem = reinterpret_cast<T*>(shdata);
    int lane = threadIdx.x % warpSize;
    int wid = threadIdx.x / warpSize;

    val = warpOp<T,op>(val, idx, size);

    if (!lane) shmem[wid] = val; // first lane of each warp fills memory

    int warps = (blockDim.x+warpSize-1)/warpSize;
    __syncthreads();
    if (threadIdx.x < warps)
        val = shmem[lane];
    // only first lane of first warp ends up with real value
    if (!wid) val = warpOp<T,op>(val, threadIdx.x, warps);

    return val;
}

template<typename T, Operation op>
__global__ // every block's sum is written to out[]
void kernel(int size, T* in, T* out) { // if called twice by reduce, in=out
    PR_assert(gridDim.x <= line/sizeof(T)); // each block writes to unique out[] index
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    T ret = reset<T,op>();
    for (int i = idx; i < size; i += blockDim.x*gridDim.x)
       ret = compare<T,op>(ret,in[i]);

    ret = blockOp<T,op>(ret, idx, size);
    if (!threadIdx.x)
        out[blockIdx.x] = ret; // block result is reused in 2nd call
}
#endif

template<typename T, Operation op = Abs>
class Reduction 
{
public:
    Reduction<T,op>() = delete;

    Reduction<T,op>(int mem = 0);

    ~Reduction<T,op>(); 

    T fetch(); // return value, waiting for kernel completion in a GPU run 

    void reduce(const T *in); // configures and calls the kernel

    bool update(T val); // compares host to val, returning true if host was updated

private:
    T *host, *in;
    int size;
#ifdef PROTO_CUDA
    T *out;
    int threads, blocks, warp;
#endif
};

template<typename T, Operation op>
Reduction<T,op>::Reduction(int mem) {
    size = mem;
#ifdef PROTO_CUDA
    protoDeviceProp prop;
    protoGetDeviceProperties(&prop, 0); // assumes system has identical GPUs
    threads = min(size,prop.maxThreadsPerBlock); // threads in a block
    blocks = (threads+size-1)/threads; // blocks in a kernel
    warp = prop.warpSize; // threads in a warp
    protoMallocHost(host,sizeof(T));
    protoMalloc(in,size*sizeof(T));
    protoMalloc(out,blocks*sizeof(T));
#else
    host = new T;
#endif
}

template<typename T, Operation op>
Reduction<T,op>::~Reduction() {
#ifdef PROTO_CUDA
    protoFreeHost(host);
    protoFree(in);
    protoFree(out);
#else
    delete host;
#endif
}

template<typename T, Operation op>
T Reduction<T,op>::fetch() {
#ifdef PROTO_CUDA
    cudaDeviceSynchronize();
    protoMemcpy(host,out,sizeof(T),protoMemcpyDeviceToHost);
#endif
    return *host;
}

template<typename T, Operation op>
void Reduction<T,op>::reduce(const T *data) {
#ifdef PROTO_CUDA
    protoMemcpy(in,data,size*sizeof(T),protoMemcpyHostToDevice);
    int shmem = (threads+warp-1)/warp; // one shmem spot per warp
    kernel<T,op><<<blocks, shmem*warp, shmem*sizeof(T)>>>(size, in, out);
    threads = min(blocks,threads); // each block in 1st kernel left a partial reduction
    shmem = (threads+warp-1)/warp;
	kernel<T,op><<<1, shmem*warp, shmem*sizeof(T)>>>(blocks, out, out); // each block made an entry in 1st call
#else
    T tmp = reset<T,op>();
    for(int it=0 ; it < size ; it++)
	    tmp = compare<T,op>(tmp,in[it]);
    *host = tmp;
#endif
}

template<typename T, Operation op>
bool Reduction<T,op>::update(T val) {
#ifdef PROTO_CUDA
    PR_error("Reduction::update shouldn't be called in CUDA run");
#else
    *host = compare<T,op>(*host, val);
#endif
    return *host == val;
}

} // end namespace Proto

#endif  // __PROTO_REDUCTION_H__
