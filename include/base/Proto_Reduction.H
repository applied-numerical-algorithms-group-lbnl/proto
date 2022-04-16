#ifndef __PROTO_REDUCTION_H__
#define __PROTO_REDUCTION_H__

//#include "Proto_gpu.H"
#include "Proto_cuda.H"
#include "Proto_macros.H"

#include <limits>
#include <typeinfo>

namespace Proto {

enum Operation { Max, Min, Abs, Sum, SumAbs };
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

template<typename T, Operation OP>
CUDA_DECORATION // every block's sum is written to out[]
T update(T last, T next) {
    T res;
    switch(OP) {
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
            break;
        case SumAbs:
            res = (next > 0 ? last + next : last - next);
            break;
    }
   return res;
}

template<typename T, Operation OP>
CUDA_DECORATION
void init(T* ptr) {
    switch (OP) {
        case Max:
            *ptr = std::numeric_limits<T>::min();
            break;
        case Min:
            *ptr = std::numeric_limits<T>::max();
            break;
        default: // Abs, Sum, SumAbs
            *ptr = static_cast<T>(0);
            break;
    }
}

template<typename T, Operation OP>
CUDA_KERNEL
void initKernel(T* ptr) 
{
    init<T, OP>(ptr);
}


#ifdef PROTO_CUDA

template<typename T, Operation OP>
__device__ // reduction within a warp, only lane 0 gets right result
T warpOp(T val, int idx, int size) {
    unsigned mask = 0xffffffff; // FULL_MASK
    if ((size/warpSize)*warpSize <= idx) // is thread in the top warp?
        mask = (1 << size%warpSize) - 1; // set bits = # of active threads in top warp
    if (idx < size) {
        for (unsigned int delta = warpSize/2; delta > 0; delta /= 2)
#ifdef PROTO_HIP
            val = update<T,OP>(val, (T)__shfl_down(val,delta)); // __shfl_down depreciated since CUDA 9.0 but __shfl_down_sync isn't available with hipcc (17/07/2020)
#else
            val = update<T,OP>(val, (T)__shfl_down_sync(mask,val,delta));
#endif
    }
    return val;
}

template<typename T, Operation OP>
__device__ // reduce within a block by storing partial warp reductions in shmem
T blockOp(T val, int idx, int size) {
    extern __shared__ __align__(sizeof(T)) unsigned char shdata[];
    T *shmem = reinterpret_cast<T*>(shdata);
    int lane = threadIdx.x % warpSize;
    int wid = threadIdx.x / warpSize;

    val = warpOp<T,OP>(val, idx, size);

    if (!lane) shmem[wid] = val; // first lane of each warp fills memory

    int warps = (blockDim.x+warpSize-1)/warpSize;
    __syncthreads();
    if (threadIdx.x < warps)
        val = shmem[lane];
    // only first lane of first warp ends up with real value
    if (!wid) val = warpOp<T,OP>(val, threadIdx.x, warps);

    return val;
}

template<typename T, Operation OP>
__global__ // every block's sum is written to out[]
void kernel(size_t size, const T* in, T* out, T* val) { // if called twice by reduce, out=nullptr
    PR_assert(gridDim.x <= line/sizeof(T)); // each block writes to unique out[] index
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    T ret = *val; // ensures that each reduction starts with result of previous one

    for (size_t i = idx; i < size; i += blockDim.x*gridDim.x)

       ret = update<T,OP>(ret,in[i]);

    ret = blockOp<T,OP>(ret, idx, size);
    if (!threadIdx.x) {
      if (gridDim.x > 1)
        out[blockIdx.x] = ret; // block result is reused in 2nd call
      else // only one block in this kernel -- always true for 2nd call
        *val = ret;
    }
}
#endif

template<typename T, Operation OP, MemType MEM>
class Reduction 
{
public:
    Reduction<T,OP,MEM>() : Reduction<T,OP,MEM>(false) {}

    // static allocation at construction, or dynamic allocation in reduce
    Reduction<T,OP,MEM>(bool dynamic);

    /// Destructor
    ~Reduction<T,OP,MEM>(); 
    
    /// Get Reduction
    /**
        Returns the computed reduction, waiting for kernel completion on GPUs
    */
    T fetch();

    /// Compute Reduction
    /**
        Calculates a reduction on the buffer <code>a_data</code> of size <code>a_size</code>
    */
    void reduce(const T *a_data, const size_t a_size); // configures and calls the kernel

    //TODO: Documentation
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

template<typename T, Operation OP, MemType MEM>
void Reduction<T,OP,MEM>::reset() {
    if (MEM == DEVICE)
    {
#ifdef PROTO_CUDA
        initKernel<T,OP><<<1,1>>>(ult);
#else
#endif
    } else {
        init<T,OP>(host);
    }
}

template<typename T, Operation OP, MemType MEM>
Reduction<T,OP,MEM>::Reduction(bool dynamic) {
    host = (T*)proto_malloc<HOST>(sizeof(T));
    if (MEM == DEVICE)
    {
#ifdef PROTO_CUDA
        // TODO: get number of streams in runtime
        //dyn = dynamic;
        dyn = false;
        protoDeviceProp prop;
        protoGetDeviceProperties(&prop, 0); // assumes system has identical GPUs
        threads = prop.maxThreadsPerBlock; // threads in a block
        warp = prop.warpSize; // threads in a warp
        ult = (T*)proto_malloc<MEM>(sizeof(T));
        //protoMallocHost(host,sizeof(T));
        //protoMalloc(MEM,ult,sizeof(T));
        initKernel<T,OP><<<1,1>>>(ult);
        if (!dyn) {
            // filling device with blocks for this kernel
            blocks = (prop.maxThreadsPerMultiProcessor/threads)*prop.multiProcessorCount; // TODO: devide by streams
            //protoMalloc(DEVICE,temp,blocks*sizeof(T)); // TODO: multiply by streams
            temp = (T*)proto_malloc<DEVICE>(blocks*sizeof(T));
        }
#else
#endif
    } else {
        init<T,OP>(host);
    }
}

template<typename T, Operation OP, MemType MEM>
Reduction<T,OP,MEM>::~Reduction() {
    proto_free<HOST>(host);
    if (MEM == DEVICE)
    {
#ifdef PROTO_CUDA
        proto_free<MEM>(ult);
        proto_free<DEVICE>(temp);
        //protoFreeHost(host);
        //protoFree(MEMTYPE_DEFAULT,ult);
        //protoFree(DEVICE,temp);
#endif
    }
}

template<typename T, Operation OP, MemType MEM>
T Reduction<T,OP,MEM>::fetch() {
    if (MEM == DEVICE)
    {
#ifdef PROTO_CUDA
        proto_memcpy<MEM, HOST>(ult, host, sizeof(T));
        pout() << "Reduction::fetch | host after memcpy: " << *host << std::endl;
        // TODO: synchronize multiple streams
        //protoMemcpy(DEVICE,host,ult,sizeof(T),protoMemcpyDeviceToHost);
        return *host;
#endif
    } else {
        return *host;
    }
}

template<typename T, Operation OP, MemType MEM>
void Reduction<T,OP,MEM>::reduce(const T *data, const size_t size) {
    if (MEM == DEVICE)
    {
#ifdef PROTO_CUDA
        if (dyn)
        {
            protoDeviceProp prop;
            protoGetDeviceProperties(&prop, 0); // assumes system has identical GPUs
            threads = min(size,(size_t)threads); // threads in a block
            blocks = (size+threads-1)/threads; // blocks in a kernel
            //protoMalloc(DEVICE,temp,blocks*sizeof(T)); // first kernel call gives one reduction per block
            temp = (T*)proto_malloc<DEVICE>(blocks*sizeof(T));
        }
        T* temp_h = (T*)proto_malloc<HOST>(blocks*sizeof(T));
        int shmem = (threads+warp-1)/warp; // warps in a block
        pout() << "Reduction::reduce | blocks: " << blocks << " | threads: " << threads << " | warp: " << warp << " | shmem: " << shmem << std::endl;
        proto_memcpy<MEM, HOST>(ult, host, sizeof(T));
        pout() << "Reduction::reduce | ult before first kernel launch: " << *host << std::endl;
        protoLaunchKernelMemAsyncGPU((kernel<T,OP>),blocks,shmem*warp,shmem*sizeof(T),
                (protoStream_t) 0, size, data, temp, ult);
        if (blocks > 1) {
            proto_memcpy<MEM, HOST>(ult, host, sizeof(T));
            pout() << "Reduction::reduce | ult before second kernel launch: " << *host << std::endl;
            proto_memcpy<MEM, HOST>(temp, temp_h, blocks*sizeof(T));
            pout() << "Reduction::reduce | temp before second kernel launch: " << std::endl;
            T temp_sum = 0;
            for (int ii = 0; ii < blocks; ii++)
            {
                pout() << temp_h[ii] << ", ";
                temp_sum += temp_h[ii];
            }
            pout() << std::endl << "Sum: " << temp_sum << std::endl;
            shmem = (blocks+warp-1)/warp; // each block in first kernel left partial reduction
            protoLaunchKernelMemAsyncGPU((kernel<T,OP>),1,shmem*warp,shmem*sizeof(T),
                    (protoStream_t) 0, blocks, temp, nullptr, ult); // updates ult
            proto_memcpy<MEM, HOST>(ult, host, sizeof(T));
            pout() << "Reduction::reduce | After both kernel launches: " << *host << std::endl;
        }
#else
        MayDay<void>::Abort("Reduction::reduce | Error: Unrecognized acceleration flag.");
#endif
    } else if (MEM == HOST)
    {
        T val = *host;
        for(int it=0 ; it < size ; it++)
            val = update<T,OP>(val,data[it]);
        *host = val;
    }
}

} // end namespace Proto

#endif  // __PROTO_REDUCTION_H__
