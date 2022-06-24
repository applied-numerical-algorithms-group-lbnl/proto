#include "Proto.H"

using namespace Proto;

void printBuffer(int* buffer, int size)
{
    for (int ii = 0; ii < size; ii++)
    {
        pout() << buffer[ii] << ", ";
    }
    pout() << std::endl;
}

int main(int argc, char** argv)
{
    #ifdef PR_MPI
    #ifdef PROTO_CUDA
    MPI_Init(&argc, &argv);
    
#if 0
    int* val_h;
    int* val_d;
    //val_h = (int*)malloc(sizeof(int));
    //cudaMalloc((void**)&val_d, sizeof(int));
    protoMalloc(HOST,   val_h, sizeof(int));
    protoMalloc(DEVICE, val_d, sizeof(int));

    *val_h = -1;
    if (procID() == 0)
    {
        *val_h = 42;
        cudaMemcpy(val_d, val_h, sizeof(int), cudaMemcpyHostToDevice);
        printf("Rank %i broadcasting %i\n", procID(), *val_h);
    } else {
        printf("Rank %i initialized with %i\n", procID(), *val_h);
    }

    MPI_Bcast(val_d, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (procID() != 0)
    {
        cudaMemcpy(val_h, val_d, sizeof(int), cudaMemcpyDeviceToHost);
        printf("Rank %i received broadcasted value %i\n", procID(), *val_h);
    }

    //free(val_h);
    //cudaFree(val_d);
    protoFree(HOST, val_h);
    protoFree(DEVICE, val_d);
    #endif
    // ALLOCATE DATA
    int* globalBuffer_h;
    int* localBuffer_h;
    int* globalBuffer_d;
    int* localBuffer_d;
    protoMalloc(HOST, globalBuffer_h, numProc()*sizeof(int));
    protoMalloc(HOST, localBuffer_h,  sizeof(int));
    protoMalloc(DEVICE, globalBuffer_d, numProc()*sizeof(int));
    protoMalloc(DEVICE, localBuffer_d,  sizeof(int));

    // INITIALIZE HOST DATA
    *localBuffer_h = (procID()+1)*7;
    for (int ii = 0; ii < numProc(); ii++)
    {
        globalBuffer_h[ii] = -1;
    }
    pout() << "proc " << procID() << " is sending " << *localBuffer_h << std::endl;
    printBuffer(globalBuffer_h, numProc());
    
    // COPY HOST -> DEVICE
    cudaMemcpy(localBuffer_d, localBuffer_h, sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(globalBuffer_d, globalBuffer_h, numProc()*sizeof(int), cudaMemcpyHostToDevice);
    pout() << "copied to device" << std::endl;
    
    // DATA MOVEMENT
    MPI_Allgather(localBuffer_d, 1, MPI_INT, globalBuffer_d, 1, MPI_INT, MPI_COMM_WORLD);
    
    // COPY DEVICE -> HOST
    cudaMemcpy(localBuffer_h, localBuffer_d, sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(globalBuffer_h, globalBuffer_d, numProc()*sizeof(int), cudaMemcpyDeviceToHost);
    pout() << "copied from device" << std::endl;
    printBuffer(globalBuffer_h, numProc());

    // FREE BUFFERS
    protoFree(HOST,     globalBuffer_h);
    protoFree(HOST,     localBuffer_h);
    protoFree(DEVICE,   globalBuffer_d);
    protoFree(DEVICE,   localBuffer_d);
    MPI_Finalize();
    #endif
    #endif
}

