
#include <mpi.h>
#include <cstdio>

int main(int argc, char* argv[])
{
  printf("Before MPI_Init\n");
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);
//  MPI_Init(&argc, &argv);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  printf("MPI rank %d\n",rank);

  MPI_Finalize();
  return 0;
}
