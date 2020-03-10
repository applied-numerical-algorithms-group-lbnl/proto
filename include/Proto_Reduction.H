#ifndef __PROTO_REDUCTION_H__
#define __PROTO_REDUCTION_H__

#include <limits>

enum Operation { Max, Min, Abs };
enum Atomic { Warp, Block, None };

constexpr size_t line = 128; // bytes in a cache line

/*template<typename T>
__device__
void add(T* ptr, T val) {
    if (sizeof(T) == sizeof(float))
        atomicAdd(static_cast<float*>(ptr), (float)val);
    else if (sizeof(T) == sizeof(double))
        atomicAdd(static_cast<double*>(ptr), (double)val);
}*/

template<typename T, Operation op>
__device__
T compare(T last, T next) {
    switch(op) {
        case Max:
            return fmax(last,next);
        case Min:
            return fmin(last,next);
        case Abs:
            return (next > 0 ? fmax(last,next) : fmax(last,-next));
        default:
            assert(0);
    }
}

template<typename T, Operation op>
__device__ // lane 0 ends up with correct value
T warpOp(T val) {
    unsigned mask = 0xffffffff;
    for (unsigned int delta = warpSize/2; delta > 0; delta /= 2)
        val = compare<T,op>(val, (T)__shfl_down_sync(mask,val,delta));
    return val;
}

template<typename T, Operation op>
__device__ 
T blockOp(T val) {
    assert(blockDim.x <= warpSize*warpSize); // no bank conflicts between warp reductions
    extern __shared__ T shmem[];
    int lane = threadIdx.x & 0x1f;
    int wid = threadIdx.x / warpSize;

    val = warpOp<T,op>(val);

    if (!lane) shmem[wid] = val; // first lane of each warp fills shmem

    __syncthreads();
    // assumes that warpSize evenly divides blockDim.x
    val = (threadIdx.x < blockDim.x/warpSize) ? shmem[lane] : 0;
    // only first lane of first warp ends up with real value
    if (!wid) val = warpOp<T,op>(val);

    return val;
}

template<typename T, Operation op, Atomic atom = Warp>
__global__ // every block's sum is written to out[]
void kernel(size_t size, T* in, T* out) { 
    assert(atom != None || gridDim.x <= line/sizeof(T)); // each block writes to unique out[] index
    T ret = out[0]; // the reset value, or the result of last kernel call
    for (int i = blockIdx.x*blockDim.x + threadIdx.x; i < size; i += blockDim.x*gridDim.x)
        ret = compare<T,op>(ret,in[i]);
    switch(atom) {
        case Warp:
            ret = warpOp<T,op>(ret);
            if (!(threadIdx.x & (warpSize-1)))
                atomicAdd(out, ret);
            break;
        case Block:
            ret = blockOp<T,op>(ret);
            if (!threadIdx.x)
                atomicAdd(out, ret);
            break;
        case None:
            ret = blockOp<T,op>(ret);
            if (!threadIdx.x)
                out[blockIdx.x] = ret;
            break;
    }
}

template<typename T, Operation op, Atomic atom = Warp>
class Reduction 
{
public:
    __host__
    Reduction<T,op,atom>() {
        cudaMallocHost(&host,sizeof(T));
        cudaMalloc(&dev,line); // 16*sizeof(Real)
        cudaEventCreate(&done);
    }

    __host__
    ~Reduction<T,op,atom>() {
        cudaFreeHost(host);
        cudaFree(dev);
        cudaEventDestroy(done);
    }

    __host__
    T fetch() {
        cudaEventSynchronize(done);
        cudaMemcpy(host,&dev[0],sizeof(T),cudaMemcpyDeviceToHost);
        return *host;
    }

    __host__
    void reset() {
        switch(op) {
            case Max:
                dev[0] = std::numeric_limits<T>::min();
            case Min:
                dev[0] = std::numeric_limits<T>::max();
            case Abs:
                dev[0] = T(0);
        }
    }

    __host__
    void reduce(T *in, const size_t size, cudaStream_t stream) {
        int warp, gpu;
        cudaGetDevice(&gpu);
        cudaDeviceGetAttribute(&warp, cudaDevAttrWarpSize, gpu);
        int blocks = std::min(line/sizeof(T), (size+warp-1)/warp);
        int threads = std::min((size_t)1024, size/blocks);
        size_t shmem = (atom==None ? threads*sizeof(T) : 0);
        kernel<T,op,atom><<<blocks,threads,shmem,stream>>>(size,in,dev);
        if (atom == None) {
            threads = min(1024,threads);
            kernel<T,op,atom><<<1,threads,0,stream>>>(blocks,dev,dev); // each block made an entry in 1st call
            cudaEventRecord(done);
        }
    }

private:
    T* host, *dev;
    cudaEvent_t done;
};

#endif // __PROTO_REDUCTION_H__
