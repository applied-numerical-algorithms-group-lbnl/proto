#ifndef __PROTO_REDUCTION_H__
#define __PROTO_REDUCTION_H__

#ifdef PROTO_CUDA
#include "implem/Proto_DisjointBoxLayout.H"
#endif

#include <limits>

namespace Proto {

enum Operation { Max, Min, Abs };
enum Atomic { Warp, Block, None };

constexpr size_t line = 128; // bytes in a cache line

template<typename T, Operation op>
CUDA_DECORATION // every block's sum is written to out[]
T compare(T last, T next) {
    switch(op) {
        case Max:
            return fmax(last,next);
        case Min:
            return fmin(last,next);
        case Abs:
            return (next > 0 ? fmax(last,next) : fmax(last,-next));
    }
}

#ifdef PROTO_CUDA

template<typename T, Operation op>
__device__ // every block's sum is written to out[]
T warpOp(T val) {
    unsigned mask = 0xffffffff;
    for (unsigned int delta = warpSize/2; delta > 0; delta /= 2)
        val = compare<T,op>(val, (T)__shfl_down_sync(mask,val,delta));
    return val;
}

template<typename T, Operation op>
__device__ // every block's sum is written to out[]
T blockOp(T val) {
    extern __shared__ T shmem[];
    int lane = threadIdx.x & 0x1f;
    int wid = threadIdx.x / warpSize;

    val = warpOp<T,op>(val);

    if (!lane) shmem[wid] = val; // first lane of each warp fills shmem

    __syncthreads();
    val = (threadIdx.x < (blockDim.x+warpSize-1)/warpSize) ? shmem[lane] : 0;
    // only first lane of first warp ends up with real value
    if (!wid) val = warpOp<T,op>(val);

    return val;
}

template<typename T, Operation op, Atomic atom = Block>
__global__ // every block's sum is written to out[]
void kernel(size_t size, T* in, T* out) { // if called twice by reduce, in=out
    PR_assert(atom != None || gridDim.x <= line/sizeof(T)); // each block writes to unique out[] index
    T ret = out[0]; // the reset value, or the result of last kernel call
    for (int i = blockIdx.x*blockDim.x + threadIdx.x; i < size; i += blockDim.x*gridDim.x)
        ret = compare<T,op>(ret,in[i]);
    switch(atom) {
        case Warp:
            ret = warpOp<T,op>(ret);
            if (!(threadIdx.x & (warpSize-1))) // first thread in a warp
                atomicAdd(out, ret);
            break;
        case Block:
            ret = blockOp<T,op>(ret);
            if (!threadIdx.x) // first thread in a block
                atomicAdd(out, ret);
            break;
        case None:
            ret = blockOp<T,op>(ret);
            if (!threadIdx.x)
                atomicAdd(&out[blockIdx.x],ret); // block result is reused in 2nd call
            break;
    }
}
#endif

template<typename T, Operation op = Abs, Atomic atom = Block>
class Reduction 
{
public:
    Reduction<T,op,atom>() : Reduction(1) {}

    Reduction<T,op,atom>(size_t mem); 

    ~Reduction<T,op,atom>(); 

    T fetch(); // return value, waiting for kernel completion in a GPU run 

    void reset(); // sets initial reduction value (depending on Operation)

    void reduce(T *in, const size_t size); // configures and calls the kernel

    bool update(T val); // compares host to val, returning true if host was updated

private:
    T *host;
#ifdef PROTO_CUDA
    T *dev;
    int thread, warp;
#endif
};

template<typename T, Operation op, Atomic atom>
Reduction<T,op,atom>::Reduction(size_t mem) {
#ifdef PROTO_CUDA
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0); // assumes system has identical GPUs
    thread = prop.maxThreadsPerBlock;
    warp = prop.warpSize;
    cudaMallocHost(&host,sizeof(T));
    if (mem != 1) {
        if (atom == None) {
            int block = prop.maxThreadsPerMultiProcessor;
            int proc = prop.multiProcessorCount;
            mem = (block/thread)*proc;
        } else {
            PR_warning("Reduction global device memory requested for an algorithm that won't use it");
            mem = 1;
        }
    } else
        cudaMalloc(&dev,mem*sizeof(T));
#else
    host = new T;
#endif
}

template<typename T, Operation op, Atomic atom>
Reduction<T,op,atom>::~Reduction() {
#ifdef PROTO_CUDA
    cudaFreeHost(host);
    cudaFree(dev);
#else
    delete host;
#endif
}

template<typename T, Operation op, Atomic atom>
T Reduction<T,op,atom>::fetch() {
#ifdef PROTO_CUDA
    cudaMemcpy(host,dev,sizeof(T),cudaMemcpyDeviceToHost);
#endif
    return *host;
}

template<typename T, Operation op, Atomic atom>
void Reduction<T,op,atom>::reset() {
    switch(op) {
        case Max:
            *host = std::numeric_limits<T>::min();
            break;
        case Min:
            *host = std::numeric_limits<T>::max();
            break;
        case Abs:
            *host = T(0);
            break;
    }
#ifdef PROTO_CUDA
    cudaMemcpy(dev,host,sizeof(T),cudaMemcpyHostToDevice);
#endif
}

template<typename T, Operation op, Atomic atom>
void Reduction<T,op,atom>::reduce(T *in, size_t size) {
#ifdef PROTO_CUDA
    int threads = min(thread,(int)size);
    int blocks = (size+threads-1)/threads;
    size_t shmem = ((threads+warp-1)/warp)*sizeof(T); // one shmem spot per warp
    cudaStream_t stream = DisjointBoxLayout::getCurrentStream();
    kernel<T,op,atom><<<blocks,threads,shmem,stream>>>(size,in,dev);
    if (atom == None) {
        threads = min(blocks,thread); // each block in 1st kernel left one partial reduction
        shmem = ((threads+warp-1)/warp)*sizeof(T);
        kernel<T,op,atom><<<1,threads,shmem,stream>>>(blocks,dev,dev); // each block made an entry in 1st call
    }
#else
    PR_error("Reduction::reduce shouldn't be called without CUDA");
#endif
}

template<typename T, Operation op, Atomic atom>
bool Reduction<T,op,atom>::update(T val) {
#ifdef PROTO_CUDA
    PR_error("Reduction::update shouldn't be called in CUDA run");
#else
    *host = compare<T,op>(*host, val);
#endif
    return *host == val;
}

} // end namespace Proto

#endif  // __PROTO_REDUCTION_H__
