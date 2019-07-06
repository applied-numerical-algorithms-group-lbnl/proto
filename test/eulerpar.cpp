#define PI 3.141592653589793
//#define DIM 3
#define NUMCELLS 64
#define NGHOST 4
//#define NUMCOMPS DIM+2

#include <algorithm>
using std::copy;
#include <vector>
using std::vector;

#include <util/Lists.hpp>
#include <util/LIKWID.hpp>

#include <mpi.h>
#include <omp.h>

#ifdef VTUNEPERF
#include <ittnotify.h>
#endif

#ifdef DATAFLOW_CODE
#include "euler_step.h"
#else
#include "Proto.H"
#include "EulerOp.H"
#endif

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

#ifndef DATAFLOW_CODE
void proto_init(double* U, Box& dbx0, BoxData<double,NUMCOMPS>& Uave, BoxData<double,NUMCOMPS>& dxdu) {
    // Setup Proto objects
#if DIM>2
    dbx0 = Box(Point(0,0,0), Point(NUMCELLS-1,NUMCELLS-1,NUMCELLS-1));
#else
    dbx0 = Box(Point(0,0), Point(NUMCELLS-1,NUMCELLS-1));
#endif
    Box dbx = dbx0.grow(NGHOST);
    Uave = BoxData<double,NUMCOMPS>(U, dbx);
    dxdu = BoxData<double,NUMCOMPS>(dbx0);
}
#endif

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
    double tsum = 0.0;
    double velmax;
    MPI_Comm comm;
    int nproc = 0;
    int pid = 0;
    int nruns = 1;
#ifdef DATAFLOW_CODE
    const char* name = "euler_step";
#else
    const char* name = "Euler::step";
    Box dbx0;
    BoxData<double,NUMCOMPS>Uave;
    BoxData<double,NUMCOMPS> dxdu;
#endif

    data_init(&U, &rhs);
#ifndef DATAFLOW_CODE
    proto_init(U, dbx0, Uave, dxdu);
#endif
    mpi_init(argc, argv, comm, &nproc, &pid);
    //nproc = omp_get_num_threads();

    if (argc > 1) {
        nruns = atoi(argv[1]);
    }

    for (unsigned i = 0; i < nruns; i++) {
    //#pragma omp parallel for private(pid)
    //for (unsigned p = 0; p < nproc; p++) {
    //pid = omp_get_thread_num();
    ptime = MPI_Wtime(); // omp_get_wtime();

#ifdef VTUNEPERF
    __SSC_MARK(0x111); // start SDE tracing, note it uses 2 underscores
    __itt_resume(); // start VTune, again use 2 underscores
#endif

#ifdef DATAFLOW_CODE
    velmax = euler_step(U, rhs);
#else
    velmax = EulerOp::step(dxdu, Uave, dbx0);
    rhs[0] = dxdu.data()[0];
#endif

#ifdef VTUNEPERF
    __itt_pause(); // stop VTune
    __SSC_MARK(0x222); // stop SDE tracing
#endif

    ptime = MPI_Wtime() - ptime;  // omp_get_wtime() - ptime;
    tsum += ptime;
    }
    
    if (pid < 1) {
        fprintf(stdout, "%s: vmax=%lf,rhs=%lf,nprocs=%d,nruns=%d,time=%lf sec)\n",
                name, velmax, rhs[0], nproc, nruns, tsum / (double) nruns);
    }

    mpi_final(comm);
    data_final(&U, &rhs);

    return 0;
}

