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

#ifdef MPI_ENABLE
#include <mpi.h>
#else
#ifdef OMP_ENABLE
#include <omp.h>
#else
#include <sys/time.h>
#endif
#endif

#ifdef VTUNE_PERF
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

void data_init(unsigned nruns, double*** U, double*** rhs) {
//    unsigned nIn = 1, nOut = 1;
//    for (unsigned d = 0; d < DIM; d++) {
//        nOut *= NUMCELLS;
//        nIn *= (NUMCELLS + 2 * NGHOST);
//    }

    unsigned nIn = 1866240;
    unsigned nOut = 1310720;

    vector<double> Uinit(nIn);
    Lists::read<double>(Uinit, DATA_FILE);

    for (unsigned i = 0; i < nruns; i++) {
        *rhs[i] = (double*) malloc(nOut * sizeof(double));
        *U[i] = (double*) malloc(nIn * sizeof(double));
        copy(Uinit.begin(), Uinit.end(), *U[i]);
    }
}

void data_final(unsigned nruns, double*** U, double*** rhs) {
    for (unsigned i = 0; i < nruns; i++) {
        free(*U[i]);
        free(*rhs[i]);
    }
}

double get_wtime() {
#ifdef MPI_ENABLE
    return MPI_Wtime();
#else
#ifdef OMP_ENABLE
    return omp_get_wtime();
#else
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double) tv.tv_sec + (((double) tv.tv_usec) * 1E-6);
#endif
#endif
}


#ifndef DATAFLOW_CODE
void proto_init(unsigned nruns, double** U, Box& dbx0, vector<BoxData<double,NUMCOMPS> >& Uave,
                vector<BoxData<double,NUMCOMPS> >& dxdu) {
    // Setup Proto objects
#if DIM>2
    dbx0 = Box(Point(0,0,0), Point(NUMCELLS-1,NUMCELLS-1,NUMCELLS-1));
#else
    dbx0 = Box(Point(0,0), Point(NUMCELLS-1,NUMCELLS-1));
#endif
    Box dbx = dbx0.grow(NGHOST);
    for (unsigned i = 0; i < nruns; i++) {
        Uave[i] = BoxData<double,NUMCOMPS>(U[i], dbx);
        dxdu[i] = BoxData<double,NUMCOMPS>(dbx0);
    }
}
#endif

#ifdef MPI_ENABLE
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
#endif

int main(int argc, char **argv) {
    double** U;
    double** rhs;
    double ptime;
    double tsum = 0.0;
    double velmax;
    int nproc = 0;
    int pid = 0;
    int nruns = 1;
#ifdef DATAFLOW_CODE
    const char* name = "euler_step";
#else
    const char* name = "Euler::step";
    Box dbx0;
    vector<BoxData<double,NUMCOMPS> > Uave;
    vector<BoxData<double,NUMCOMPS> > dxdu;
#endif

#ifdef MPI_ENABLE
    MPI_Comm comm;
    mpi_init(argc, argv, comm, &nproc, &pid);
#else
#ifdef OMP_ENABLE
    nproc = omp_get_num_threads();
#else
    nproc = 1;
#endif
#endif

    if (argc > 1) {
        nruns = atoi(argv[1]);
    }

    data_init(nruns, &U, &rhs);
#ifndef DATAFLOW_CODE
    Uave.resize(nruns);
    dxdu.resize(nruns);
    proto_init(nruns, U, dbx0, Uave, dxdu);
#endif

#ifdef LIKWID_PERF
    LIKWID_MARKER_INIT
#endif

    for (unsigned i = 0; i < nruns; i++) {
#ifdef OMP_ENABLE
    #pragma omp parallel for private(pid)
    for (unsigned p = 0; p < nproc; p++) {
    pid = omp_get_thread_num();
#endif

    ptime = get_wtime();

#ifdef LIKWID_PERF
    LIKWID_MARKER_START(name);
#endif

#ifdef VTUNE_PERF
    __SSC_MARK(0x111); // start SDE tracing, note it uses 2 underscores
    __itt_resume(); // start VTune, again use 2 underscores
#endif

#ifdef DATAFLOW_CODE
    velmax = euler_step(U[i], rhs[i]);
#else
    velmax = EulerOp::step(dxdu[i], Uave[i], dbx0);
    *rhs[0] = dxdu[i].data()[0];
#endif

#ifdef LIKWID_PERF
    LIKWID_MARKER_STOP(name);
#endif

#ifdef VTUNE_PERF
    __itt_pause(); // stop VTune
    __SSC_MARK(0x222); // stop SDE tracing
#endif

    ptime = get_wtime() - ptime;

#ifdef OMP_ENABLE
    }
#endif
    tsum += ptime;
    }

#ifdef LIKWID_PERF
    LIKWID_MARKER_CLOSE
#endif

    if (pid < 1) {
        fprintf(stdout, "%s: vmax=%lf,rhs=%lf,nprocs=%d,nruns=%d,time=%lf\n",
                name, velmax, *rhs[0], nproc, nruns, tsum / (double) nruns);
    }

    data_final(nruns, &U, &rhs);
#ifdef MPI_ENABLE
    mpi_final(comm);
#endif

    return 0;
}

