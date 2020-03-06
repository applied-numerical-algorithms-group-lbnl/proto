#include <limits>

enum Operation { Max, Min, Abs };
enum Atomic { Warp, Block, None };

constexpr size_t line = 128; // bytes in a cache line

template<typename T, Operation op, Atomic atom = Warp>
class Reduction 
{
public:
    Reduction<T,op>() : {
        cudaMallocHost(&host,sizeof(T));
        cudaMalloc(&dev,const_cast<size_t>(line));
        cudaMemcpyToSymbol(in
        cudaEventCreate(&done);
    }

    ~Reduction<T>() : {
        cudaFreeHost(host);
        cudaFree(dev);
        cudaEventDestroy(done);
    }

    T fetch() {
        cudaEventSynchronize(done);
        cudaMemcpy(host,&dev[0],sizeof(T),cudaMemcpyDeviceToHost);
        return *host;
    }

    void reset() {
        switch(op) {
            case max:
                dev[0] = std::numeric_limits<T>::min();
            case min:
                dev[0] = std::numeric_limits<T>::max();
            case abs:
                dev[0] = T(0);
        }
    }

    T* ptr() {
        return dev;
    }

    __host__
    void reduce(int *in, const size_t size, cudaStream_t stream) {
        int blocks = min(line/sizeof(T), (size+warpSize-1)/warpSize);
        int threads = min(1024, size/blocks);
        kernel<<<blocks,threads,stream>>>(size,in);

        threads = min(1024,threads);
        kernel<<<1,threads, stream>>>(blocks); // each block made an entry in last call
        cudaEventRecord(done);
    }

private:

    __device__ // every block's sum is written to out[]
    void kernel(size_t size, T* ptr = dev) { // use default argument 
        assert(atom != none || gridDim.x <= line/sizeof(T)); // each block writes to unique out[] index
        ret = dev[0]; // the reset value, or the result of last kernel call
        for (int i = blockIdx.x*blockDim.x + threadIdx.x; i < size; i += blockDim.x*gridDim.x)
            ret = compare(ret,ptr[i]);
        switch(atom) {
            case warp:
                ret = warpOp(ret);
                if (!(threadIdx.x & (warpSize-1)))
                    atomicAdd(dev, sum);
                break;
            case block:
                ret = blockOp(ret);
                if (!threadIdx.x)
                    atomicAdd(dev, sum);
                break;
            case none:
                ret = blockOp(ret);
                if (!threadIdx.x)
                    dev[blockIdx.x] = ret;
                break;
        }
    }
    
    __device__
    T fn(T last, T next) {
        switch(op) {
            case max:
                return max(last,next);
            case min:
                return min(last,next);
            case abs:
                return (next > 0 ? max(last,next) : max(last,-next);
            default:
                assert(0);
        }
    }

    __device__ // lane 0 ends up with correct value
    T warpOp(T val) {
        unsigned mask = 0xffffffff;
        for (unsigned int delta = warpSize/2; delta > 0; delta /= 2)
            val = fn(val, __shfl_down_sync(mask,val,delta));
        return val;
    }

    __device__ 
    T blockOp(T val) {
        assert(blockDim.x <= warpSize*warpSize); // no bank conflicts between warp reductions
        static __shared__ int shared[blockDim.x];
        int lane = threadIdx.x & 0x1f;
        int wid = threadIdx.x / warpSize;

        val = warpOp(val);

        if (!lane) shared[wid] = val; // first lane of each warp fills shmem

        __syncthreads();
        // assumes that warpSize evenly divides blockDim.x
        val = (threadIdx.x < blockDim.x/warpSize) ? shared[lane] : 0;
        // only first lane of first warp ends up with real value
        if (!wid) val = warpOp(val);

        return val;
    }

    T* host, *dev;
    cudaEvent_t done;
};
