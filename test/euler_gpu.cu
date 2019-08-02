#define PI 3.141592653589793
#define DIM 3
#define NUMCELLS 64
#define NGHOST 4
//#define NUMCOMPS DIM+2

#include <algorithm>
using std::copy;
#include <vector>
using std::vector;

#include <Lists.hpp>
#include <LIKWID.hpp>
#include <cuda_funcs.h>

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

int main(int argc, char **argv) {
    double* U;
    double* rhs;
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
//    cuda_print(cuda, stdout);

    data_init(&U, &rhs);
#ifndef DATAFLOW_CODE
    proto_init(U, dbx0, Uave, dxdu);
#endif

    if (argc > 1) {
        nruns = atoi(argv[1]);
    }

    for (unsigned i = 0; i < nruns; i++) {
    cuda_profile_start(cuda);

#ifdef DATAFLOW_CODE
    velmax = euler_step(U, rhs);
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

    data_final(&U, &rhs);
    cuda_del(cuda);

    return 0;
}

