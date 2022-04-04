#include "Proto.H"

using namespace Proto;

__global__ void f_foo(int* a)
{
    int id = blockIdx.x*blockDim.x + threadIdx.x;
    a[id] = 7;
}

int main(int argc, char** argv)
{
    #ifdef PR_MPI
    MPI_Init(&argc, &argv);
    #endif

    size_t size = 8*sizeof(int);

    int* hostBuffer_0 = (int*)std::malloc(size);
    int* hostBuffer = (int*)proto_malloc<HOST>(size);
    MemType mem;
    std::cout << "Checking manual host buffer: " << hostBuffer_0 << std::endl;
    for (int ii = 0; ii < 8; ii++) { hostBuffer_0[ii] = 0; }
    mem = pointerMemType(hostBuffer_0);
    std::cout << "buffer memtype: " << mem << std::endl; 
    std::cout << "success." << std::endl;
    std::cout << "Checking proto host buffer: " << hostBuffer << std::endl;
    for (int ii = 0; ii < 8; ii++) { hostBuffer[ii] = 0; }
    mem = pointerMemType(hostBuffer);
    std::cout << "buffer memtype: " << mem << std::endl; 
    std::cout << "success." << std::endl;

    proto_free<HOST>(hostBuffer);

    #ifdef PROTO_CUDA
    int* deviBuffer_0;
    cudaMalloc((void**)&deviBuffer_0, size);
    int* deviBuffer = (int*)proto_malloc<DEVICE>(size);

    std::cout << "Checking manual device buffer" << std::endl;
    f_foo<<<1,8>>>(deviBuffer_0);
    mem = pointerMemType(deviBuffer_0);
    std::cout << "buffer memtype: " << mem << std::endl; 
    std::cout << "success." << std::endl;
    std::cout << "Checking proto device buffer" << std::endl;
    f_foo<<<1,8>>>(deviBuffer);
    mem = pointerMemType(deviBuffer);
    std::cout << "buffer memtype: " << mem << std::endl; 
    std::cout << "success." << std::endl;

    proto_free<DEVICE>(deviBuffer);
    #endif

    #ifdef PR_MPI
    MPI_Finalize();
    #endif
}
