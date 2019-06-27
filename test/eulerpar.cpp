#define PI 3.141592653589793
#define DIM 3
#define NUMCELLS 64
#define NGHOST 4
#define NUMCOMPS DIM+2

#include <algorithm>
using std::copy;
#include <vector>
using std::vector;

#include <util/Lists.hpp>
#include <util/LIKWID.hpp>

#include <mpi.h>
//#include <omp.h>

#ifdef VTUNEPERF
#include <ittnotify.h>
#endif

#include "euler_step.h"
#if DIM>2
#define DATA_FILE "data/Uin_3d.csv"
#else
#define DATA_FILE "data/Uin_2d.csv"
#endif

void data_init(double** U, double** rhs) {
    unsigned nIn = 1, nOut = 1;
    for (unsigned d = 0; d < DIM; d++) {
        nOut *= NUMCELLS;
        nIn *= (NUMCELLS + 2 * NGHOST);
    }

    unsigned rhs_size = 1310720;
    *rhs = (double*) malloc(rhs_size * sizeof(double));
    unsigned Usize = 1866240;
    *U = (double*) malloc(Usize * sizeof(double));
    vector<double> Uinit(Usize);
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
    double ptime;
    MPI_Comm comm;
    int nproc = 0;
    int pid = 0;

    data_init(&U, &rhs);
    mpi_init(argc, argv, comm, &nproc, &pid);
    //nproc = omp_get_num_threads();

    //#pragma omp parallel for private(pid)
    //for (unsigned p = 0; p < nproc; p++) {
    //pid = omp_get_thread_num();
    if (pid < 1) {
        ptime = MPI_Wtime(); // omp_get_wtime();
    }

#ifdef VTUNEPERF
    __SSC_MARK(0x111); // start SDE tracing, note it uses 2 underscores
    __itt_resume(); // start VTune, again use 2 underscores
#endif
    double velmax = euler_step(U, rhs);
#ifdef VTUNEPERF
    __itt_pause(); // stop VTune
    __SSC_MARK(0x222); // stop SDE tracing
#endif
    if (pid < 1) {
        ptime = MPI_Wtime() - ptime;  // omp_get_wtime() - ptime;
        fprintf(stderr, "euler_step: vmax=%lf (%lf sec)\n", velmax, ptime);
    }
    //}

    mpi_final(comm);
    data_final(&U, &rhs);

    return 0;
}

