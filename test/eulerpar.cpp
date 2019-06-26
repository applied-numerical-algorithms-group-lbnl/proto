#define PI 3.141592653589793
#define NGHOST 4
#define NCOMP DIM+2

#include <algorithm>
using std::copy;
#include <vector>
using std::vector;
#include <util/Lists.hpp>

#include <mpi.h>

#include "euler_step.h"
#if DIM>2
#define DATA_FILE "data/Uin_3d.csv"
#else
#define DATA_FILE "data/Uin_2d.csv"
#endif

void data_init(double** U, double** rhs) {
    unsigned nIn = 1, nOut = 1;
    for (unsigned d = 0; d < DIM; d++) {
        nIn *= NUMCELLS + 2 * NGHOST;
        nOut *= NUMCELLS;
    }

    *rhs = (double*) malloc(nOut * NCOMP * sizeof(double));
    *U = (double*) malloc(nIn * NCOMP * sizeof(double));

    vector<double> Uinit(nIn * NCOMP);
    Lists::read<double>(Uinit, DATA_FILE);
    copy(Uinit.begin(), Uinit.end(), *U);
}

void data_final(double** U, double** rhs) {
    free(*rhs);
    free(*U);
}

void mpi_init(int argc, char **argv, MPI_Comm& comm, int* nproc, int* pid) {
    MPI_Init(&argc, &argv);
    MPI_Comm_dup(MPI_COMM_WORLD, &comm);
    MPI_Comm_size(comm, nproc);
    MPI_Comm_rank(comm, pid);
}

void mpi_final(MPI_Comm& comm) {
    MPI_Comm_free(&comm);
    MPI_Finalize();
}

int main(int argc, char **argv) {
    double* U;
    double* rhs;
    MPI_Comm comm;
    int nproc = 0;
    int pid = 0;

    data_init(&U, &rhs);
    mpi_init(argc, argv, comm, &nproc, &pid);

    double pstart = MPI_Wtime();
    double velmax = euler_step(U, rhs);
    fprintf(stderr, "euler_step (%d/%d): %lf sec\n", pid, nproc, MPI_Wtime() - pstart);

    mpi_final(comm);
    data_final(&U, &rhs);

    return 0;
}
