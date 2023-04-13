#include "Proto.H"

using namespace Proto;

__global__ void f_initSrc(int* a)
{
    int id = blockIdx.x*blockDim.x + threadIdx.x;
    a[id] = (id+1)*2;
}

__global__ void f_initDst(int* a)
{
    int id = blockIdx.x*blockDim.x + threadIdx.x;
    a[id] = (id+1)*3;
}

int main(int argc, char** argv)
{
    #ifdef PR_MPI
    MPI_Init(&argc, &argv);
    #endif

    bool PASS = true;

    int size = 8;
    size_t numBytes = size*sizeof(int);
    
    int* srcBuffer = (int*)proto_malloc<HOST>(numBytes);
    int* dstBuffer = (int*)proto_malloc<HOST>(numBytes);

    for (int ii = 0; ii < size; ii++)
    {
        srcBuffer[ii] = (ii+1)*3;
        dstBuffer[ii] = -1;
    }
    proto_memcpy<HOST, HOST>(srcBuffer, dstBuffer, numBytes);

    for (int ii = 0; ii < size; ii++)
    {
        PASS &= (srcBuffer[ii] == dstBuffer[ii]);
        if (!PASS)
        {
            std::cout << "HOST -> HOST test failed. src = " << srcBuffer[ii];
            std::cout << " | dst = " << dstBuffer[ii] << std::endl;
        }
    }

    #ifdef PROTO_ACCEL

    int* srcBuffer_h = (int*)proto_malloc<HOST>(numBytes);
    int* srcBuffer_d = (int*)proto_malloc<DEVICE>(numBytes);
    int* dstBuffer_d = (int*)proto_malloc<DEVICE>(numBytes);
    int* dstBuffer_h = (int*)proto_malloc<HOST>(numBytes);

    for (int ii = 0; ii < size; ii++)
    {
        srcBuffer_h[ii] = (ii+1);
        dstBuffer_h[ii] = -1;
    }
    f_initSrc<<<1,size>>>(srcBuffer_d);
    f_initDst<<<1,size>>>(dstBuffer_d);

    proto_memcpy<HOST, DEVICE>(srcBuffer_h, srcBuffer_d, numBytes);
    cudaDeviceSynchronize();
    proto_memcpy<DEVICE, DEVICE>(srcBuffer_d, dstBuffer_d, numBytes);
    cudaDeviceSynchronize();
    proto_memcpy<DEVICE, HOST>(dstBuffer_d, dstBuffer_h, numBytes);
    cudaDeviceSynchronize();
    
    for (int ii = 0; ii < size; ii++)
    {
        PASS &= (srcBuffer_h[ii] == dstBuffer_h[ii]);
        if (!PASS)
        {
            std::cout << "HOST <--> DEVICE test failed. src = " << srcBuffer_h[ii];
            std::cout << " | dst = " << dstBuffer_h[ii] << std::endl;
        }
    }

    #endif

    if (PASS)
    {
        std::cout << "ALL TESTS PASSED" << std::endl;
    }

    #ifdef PR_MPI
    MPI_Finalize();
    #endif
}
