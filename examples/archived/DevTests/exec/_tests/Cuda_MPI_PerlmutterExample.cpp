#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <cuda_runtime.h>

int main(int argc, char *argv[]) {
    int myrank;
    float *val_device, *val_host;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    val_host = (float*)malloc(sizeof(float));
    cudaMalloc((void **)&val_device, sizeof(float));

    *val_host = -1.0;
    if (myrank != 0) {
      printf("%s %d %s %f\n", "I am rank", myrank, "and my initial value is:", *val_host);
    }

    if (myrank == 0) {
        *val_host = 42.0;
        cudaMemcpy(val_device, val_host, sizeof(float), cudaMemcpyHostToDevice);
        printf("%s %d %s %f\n", "I am rank", myrank, "and will broadcast value:", *val_host);
    }

    MPI_Bcast(val_device, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

    if (myrank != 0) {
      cudaMemcpy(val_host, val_device, sizeof(float), cudaMemcpyDeviceToHost);
      printf("%s %d %s %f\n", "I am rank", myrank, "and received broadcasted value:", *val_host);
    }

    cudaFree(val_device);
    free(val_host);

    MPI_Finalize();
    return 0;
}
