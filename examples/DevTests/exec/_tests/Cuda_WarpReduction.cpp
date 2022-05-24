#include "Proto.H"
#include "InputParser.H"

#define FULL_MASK 0xffffffff;

using namespace Proto;

template<typename T>
CUDA_DECORATION
T f_update(T prev, T next)
{
    return prev + next;
}

template<typename T>
__device__
void f_warpOp(T val, int idx, size_t size)
{
    unsigned mask = FULL_MASK;
    if ((size / warpSize) * warpSize <= idx)
    {
        mask  = (1 << size%warpSize) - 1;
    }
    for (unsigned int delta = warpSize / 2; delta > 0; delta /= 0)
    {
        val = f_update<T>(val, (T)__shfl_down_sync(mask, val, delta));
    }
    return val;
}

template<typename T>
__device__
void f_blockOp(T val, int idx, size_t size)
{
    extern __shared__ __align__(sizeof(T)) unsigned char shdata[];
    T* shmem = reinterpret_cast<T*>(shdata);
    int lane = threadIdx.x % warpSize;
    int warpID = threadIdx.x / warpSize;

    val = f_warpOp<T>(val, idx, size);

    if (lane == 0)
    {
        shmem[warpID] = val;
    }

}

template<typename T>
__global__
void f_reductionKernel(size_t size, const T* input, T* output, T* value)
{
    int idx = blockIdx.x*blockDim.x + threadIdx.x; // get thread id
    T ret = *value;
    for (size_t ii = idx; ii < size; ii += blockDim.x*gridDim.x)
    {
        ret = f_update<T>(ret, input[ii]);
    }
    if (threadIdx.x == 0)
    {
        if (gridDim.x > 1)
        {
            output[blockIdx.x] = ret;
        } else {
            *value = ret;
        }
    }
}

int main(int argc, char** argv)
{
#ifdef PR_MPI
    MPI_Init(&argc, &argv);
#ifdef PROTO_CUDA
    HDF5Handler h5;

    int domainSize = 64;

    InputArgs args;
    args.add("domainSize",      domainSize);
    args.parse(argc, argv);
    args.print();

    size_t numBytes = domainSize*sizeof(int);

    protoDeviceProp prop;
    protoGetDeviceProperties(&prop, 0);
    int threadsPerBlock = prop.maxThreadsPerBlock;
    int threadsPerWarp = prop.warpSize;
    int warpsPerBlock = (threadsPerBlock + threadsPerWarp - 1) / threadsPerWarp;
    int numBlocks = (prop.maxThreadsPerMultiProcessor / threadsPerBlock)*prop.multiProcessorCount;
    
    std::cout << "blocks: " << numBlocks;
    std::cout << " | threadsPerWarp: " << threadsPerWarp;
    std::cout << " | threadsPerBlock: " << threadsPerBlock;
    std::cout << " | warpsPerBlock: " << warpsPerBlock << std::endl;

    int* data_d =  (int*)proto_malloc<DEVICE>(numBytes);
    int* data_h =  (int*)proto_malloc<HOST>(numBytes);
    int* temp_d =  (int*)proto_malloc<DEVICE>(numBlocks*sizeof(int));
    int* temp_h =  (int*)proto_malloc<HOST>(numBlocks*sizeof(int));
    int* value_d = (int*)proto_malloc<DEVICE>(sizeof(int));
    int* value_h = (int*)proto_malloc<HOST>(sizeof(int));

    for (int ii = 0; ii < domainSize; ii++)
    {
        data_h[ii] = 1;
    }

    proto_memcpy<HOST, DEVICE>(data_h, data_d, numBytes);
    
    f_reductionKernel<int><<<numBlocks, threadsPerBlock>>>(domainSize, data_d, temp_d, value_d);
    
    proto_memcpy<DEVICE, HOST>(temp_d,  temp_h,  numBlocks*sizeof(int));
    proto_memcpy<DEVICE, HOST>(value_d, value_h, sizeof(int));

    std::cout << "Temp Data" << std::endl;
    for (int ii = 0; ii < numBlocks; ii++)
    {
        std::cout << temp_h[ii] << ", ";
    }
    std::cout << std::endl;
    std::cout << "Value Data" << std::endl;
    std::cout << *value_h << std::endl;
#endif
    MPI_Finalize();
#endif
}

