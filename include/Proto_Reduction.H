#ifndef __PROTO_REDUCTION_H__
#define __PROTO_REDUCTION_H__

#include "Proto_gpu.H"
#include "Proto_cuda.H"


#include <limits>
#include <typeinfo>

namespace Proto {

enum Operation { Max, Min, Abs, Sum };
enum Atomic { Warp, Block, None };

constexpr int line = 128; // bytes in a cache line

#ifdef PROTO_CUDA // namespace collision in host builds
template<typename T>
CUDA_DECORATION
static constexpr T max(T a, T b) { return a > b ? a : b; }
template<typename T>
CUDA_DECORATION
static constexpr T min(T a, T b) { return a < b ? a : b; }
#endif

template<typename T, Operation op>
CUDA_DECORATION // every block's sum is written to out[]
T update(T last, T next) {
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
        case Sum:
            res = last + next;
    }
   return res;
}

template<typename T, Operation op>
CUDA_KERNEL
void init(T* ptr) {
    switch (op) {
        case Max:
            *ptr = std::numeric_limits<T>::min();
            break;
        case Min:
            *ptr = std::numeric_limits<T>::max();
            break;
        default: // Abs and Sum
            *ptr = static_cast<T>(0);
            break;
    }
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
        val = update<T,op>(val, (T)__shfl_down(val,delta)); // __shfl_down depreciated since CUDA 9.0 but __shfl_down_sync isn't available with hipcc (17/07/2020)
#else
        val = update<T,op>(val, (T)__shfl_down_sync(mask,val,delta));
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
void kernel(size_t size, const T* in, T* out, T* val) { // if called twice by reduce, out=nullptr
    PR_assert(gridDim.x <= line/sizeof(T)); // each block writes to unique out[] index
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    T ret = *val; // ensures that each reduction starts with result of previous one

    for (size_t i = idx; i < size; i += blockDim.x*gridDim.x)

       ret = update<T,op>(ret,in[i]);

    ret = blockOp<T,op>(ret, idx, size);
    if (!threadIdx.x) {
      if (gridDim.x > 1)
        out[blockIdx.x] = ret; // block result is reused in 2nd call
      else // only one block in this kernel -- always true for 2nd call
        *val = ret;
    }
}
#endif

template<typename T, Operation op = Abs>
class Reduction 
{
public:
    Reduction<T,op>() : Reduction<T,op>(false) {}

    // static allocation at construction, or dynamic allocation in reduce
    Reduction<T,op>(bool dynamic);

    ~Reduction<T,op>(); 

    T fetch(); // return value, waiting for kernel completion in a GPU run 

    void reduce(const T *data, const size_t size); // configures and calls the kernel

    void reset(); // initializes pointer based on operation

private:
    T *host;
    int size;
#ifdef PROTO_CUDA
    bool dyn;
    T *temp, *ult;
    int threads, blocks, warp;
#endif
};

template<typename T, Operation op>
void Reduction<T,op>::reset() {
#ifdef PROTO_CUDA
    init<T,op><<<1,1>>>(ult);
#else
    init<T,op>(host);
#endif
}

template<typename T, Operation op>
Reduction<T,op>::Reduction(bool dynamic) {
#ifdef PROTO_CUDA
// TODO: get number of streams in runtime
    dyn = dynamic;
    protoDeviceProp prop;
    protoGetDeviceProperties(&prop, 0); // assumes system has identical GPUs
    threads = prop.maxThreadsPerBlock; // threads in a block
    warp = prop.warpSize; // threads in a warp
    protoMallocHost(host,sizeof(T));
    protoMalloc(ult,sizeof(T));
    init<T,op><<<1,1>>>(ult);
    if (!dyn) {
        // filling device with blocks for this kernel
        blocks = (prop.maxThreadsPerMultiProcessor/threads)*prop.multiProcessorCount; // TODO: devide by streams
        protoMalloc(temp,blocks*sizeof(T)); // TODO: multiply by streams
    }
#else
    host = new T;
    init<T,op>(host);
#endif
}

template<typename T, Operation op>
Reduction<T,op>::~Reduction() {
#ifdef PROTO_CUDA
    protoFreeHost(host);
    protoFree(ult);
    protoFree(temp);
#else
    delete host;
#endif
}

template<typename T, Operation op>
T Reduction<T,op>::fetch() {
#ifdef PROTO_CUDA
    // TODO: synchronize multiple streams
    protoMemcpy(host,ult,sizeof(T),protoMemcpyDeviceToHost);
#endif
    return *host;
}

template<typename T, Operation op>
void Reduction<T,op>::reduce(const T *data, const size_t size) {
#ifdef PROTO_CUDA
    if (dyn) {
        protoDeviceProp prop;
        protoGetDeviceProperties(&prop, 0); // assumes system has identical GPUs
        threads = min(size,(size_t)threads); // threads in a block
        blocks = (size+threads-1)/threads; // blocks in a kernel
        protoMalloc(temp,blocks*sizeof(T)); // first kernel call gives one reduction per block
    }
    int shmem = (threads+warp-1)/warp; // warps in a block

    kernel<T,op><<<blocks,shmem*warp,shmem*sizeof(T)>>>(size, data, temp, ult);
    if (blocks > 1) {
        shmem = (blocks+warp-1)/warp; // each block in first kernel left partial reduction
        kernel<T,op><<<1,shmem*warp,shmem*sizeof(T)>>>(blocks, temp, nullptr, ult); // updates ult
    }

#else
    T val = *host;
    for(int it=0 ; it < size ; it++)
	    val = update<T,op>(val,data[it]);
    *host = val;
#endif
}

} // end namespace Proto

#endif  // __PROTO_REDUCTION_H__
