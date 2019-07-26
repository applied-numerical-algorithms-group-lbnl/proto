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

#include <util/cuda_funcs.h>

#ifdef DATAFLOW_CODE
#include "euler_step_3d_dev_gpu.h"
#else
#include "Proto.H"
#include "EulerOp.H"
#endif

#if DIM>2
#define DATA_FILE "data/Uin_3d.csv"
#else
#define DATA_FILE "data/Uin_2d.csv"
#endif

void data_init(double** h_U, double** h_rhs, double** d_U, double** d_rhs) {
    unsigned nIn = 1, nOut = 1;
    for (unsigned d = 0; d < DIM; d++) {
        nOut *= NUMCELLS;
        nIn *= (NUMCELLS + 2 * NGHOST);
    }

    unsigned rhs_size = 1310720;
    *h_rhs = (double*) malloc(rhs_size * sizeof(double));
    cuda_malloc((void**) d_rhs, rhs_size);

    unsigned Usize = 1866240;
    *h_U = (double*) malloc(Usize * sizeof(double));
    cuda_malloc((void**) d_U, Usize);

    vector<double> Uinit(Usize);
    Lists::read<double>(Uinit, DATA_FILE);
    copy(Uinit.begin(), Uinit.end(), *h_U);
    cuda_copy_device(*h_U, *d_U, Usize);
}

void data_final(double** h_U, double** h_rhs, double** d_U, double** d_rhs) {
    free(*h_rhs);
    free(*h_U);
    cuda_free(*d_rhs);
    cuda_free(*d_U);
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

int main(int argc, char **argv) {
    double* U;
    double* rhs;
    double* d_U;
    double* d_rhs;
    double ptime;
    double tsum = 0.0;
    double velmax;
    int nproc = 1;
    int pid = 0;
    int nruns = 1;

#ifdef DATAFLOW_CODE
    const char* name = "euler_step";
#else
    const char* name = "Euler::step";
    Box dbx0;
    BoxData<double,NUMCOMPS> Uave;
    BoxData<double,NUMCOMPS> dxdu;
#endif

    cuda_t* cuda = cuda_new();

    data_init(&U, &rhs, &d_U, &d_rhs);
#ifndef DATAFLOW_CODE
    proto_init(U, dbx0, Uave, dxdu);
#endif

    if (argc > 1) {
        nruns = atoi(argv[1]);
    }

    for (unsigned i = 0; i < nruns; i++) {
    cuda_profile_start(&start, &stop);

#ifdef DATAFLOW_CODE
    velmax = euler_step(d_U, d_rhs);
#else
    velmax = EulerOp::step(dxdu, Uave, dbx0);
    rhs[0] = dxdu.data()[0];
#endif

    ptime = cuda_profile_stop(cuda);

    tsum += ptime;
    }

    if (pid < 1) {
        fprintf(stdout, "%s: vmax=%lf,rhs=%lf,nprocs=%d,nruns=%d,time=%lf\n",
                name, velmax, rhs[0], nproc, nruns, tsum / (double) nruns);
    }

    data_final(&U, &rhs, &d_U, &d_rhs);
    cuda_del(cuda);

    return 0;
}

