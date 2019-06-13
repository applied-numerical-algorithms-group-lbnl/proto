#include <cstring>
#include <string>
using std::string;
using std::to_string;

#include <gtest/gtest.h>
using namespace testing;
#include <BenchmarkTest.hpp>
//using test::BenchmarkTest;

#define PI 3.141592653589793
#define NGHOST 4
#define NUMCELLS 16

#include "Proto.H"
#include "EulerRK4.H"
#include "EulerOp.H"
#include "Proto_DataFlow.H"

#include "euler_step.h"

namespace test {
    class EulerTest : public BenchmarkTest {
    protected:
        EulerTest() {
            _name = "EulerTest";
        }

        virtual void SetUp(initializer_list<string> args) {
            _tstop= 1.0;
            _maxstep = 10;
            _outputInterval = -1;
            _nx = NUMCELLS;

            vector<string> argv(args.begin(), args.end());
            unsigned argc = argv.size();

            for(int iarg = 0; iarg < argc-1; iarg++) {
                if(strcmp(argv[iarg].c_str(),"-n") == 0) {
                    _nx = atoi(argv[iarg+1].c_str());
                } else if(strcmp(argv[iarg].c_str(), "-m") == 0) {
                    _maxstep = atoi(argv[iarg+1].c_str());
                } else if(strcmp(argv[iarg].c_str(), "-o") == 0) {
                    _outputInterval = atoi(argv[iarg+1].c_str());
                } else if(strcmp(argv[iarg].c_str(),"-t") == 0) {
                    _tstop = atof(argv[iarg+1].c_str());
                }
            }

            _nOut = NUMCELLS;
            _nIn = NUMCELLS + 2 * NGHOST;
            _nVec = NUMCELLS + NGHOST;

//            _nAve = (NUMCELLS + NGHOST/2) * (NUMCELLS + NGHOST + NGHOST/2);
//            _nAve2 = (NUMCELLS + 1) * (NUMCELLS + NGHOST + NGHOST/2);
//            _nFlux = (NUMCELLS + 1) * (NUMCELLS + NGHOST);

            for (unsigned d = 1; d < DIM; d++) {
                _nOut *= _nOut;
                _nIn *= _nIn;
                _nVec *= _nVec;
            }

            _rhs = (double*) calloc(_nOut * NUMCOMPS, sizeof(double));
            _Uin = (double*) calloc(_nIn * NUMCOMPS, sizeof(double));
            _rhs_out = (double*) calloc(_nOut * NUMCOMPS, sizeof(double));

            _retval = 2.76852;

            double Uinit[] = {1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,1.29546,1.14763,1.14763,1.29546,1.50454,1.65237,1.65237,1.50454,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,0.458016,0.405748,0.405748,0.458016,0.531934,0.584202,0.584202,0.531934,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941,2.40059,2.01252,2.01252,2.40059,2.94941,3.33748,3.33748,2.94941};
            for (unsigned i = 0; i < _nIn * NUMCOMPS; i++) {
                _Uin[i] = Uinit[i];
            }

            double rhs[] = {-0.0703309,-0.0298238,0.028145,0.0696183,0.0703166,0.0298381,-0.0281307,-0.0696326,-0.0703309,-0.0298238,0.028145,0.0696183,0.0703166,0.0298381,-0.0281307,-0.0696326,-0.0703309,-0.0298238,0.028145,0.0696183,0.0703166,0.0298381,-0.0281307,-0.0696326,-0.0703309,-0.0298238,0.028145,0.0696183,0.0703166,0.0298381,-0.0281307,-0.0696326,-0.0703309,-0.0298238,0.028145,0.0696183,0.0703166,0.0298381,-0.0281307,-0.0696326,-0.0703309,-0.0298238,0.028145,0.0696183,0.0703166,0.0298381,-0.0281307,-0.0696326,-0.0703309,-0.0298238,0.028145,0.0696183,0.0703166,0.0298381,-0.0281307,-0.0696326,-0.0703309,-0.0298238,0.028145,0.0696183,0.0703166,0.0298381,-0.0281307,-0.0696326,-0.0703309,-0.0298238,0.028145,0.0696183,0.0703166,0.0298381,-0.0281307,-0.0696326,-0.0703309,-0.0298238,0.028145,0.0696183,0.0703166,0.0298381,-0.0281307,-0.0696326,-0.0703309,-0.0298238,0.028145,0.0696183,0.0703166,0.0298381,-0.0281307,-0.0696326,-0.0703309,-0.0298238,0.028145,0.0696183,0.0703166,0.0298381,-0.0281307,-0.0696326,-0.0703309,-0.0298238,0.028145,0.0696183,0.0703166,0.0298381,-0.0281307,-0.0696326,-0.0703309,-0.0298238,0.028145,0.0696183,0.0703166,0.0298381,-0.0281307,-0.0696326,-0.0703309,-0.0298238,0.028145,0.0696183,0.0703166,0.0298381,-0.0281307,-0.0696326,-0.0703309,-0.0298238,0.028145,0.0696183,0.0703166,0.0298381,-0.0281307,-0.0696326,-0.0703309,-0.0298238,0.028145,0.0696183,0.0703166,0.0298381,-0.0281307,-0.0696326,-0.0703309,-0.0298238,0.028145,0.0696183,0.0703166,0.0298381,-0.0281307,-0.0696326,-0.0703309,-0.0298238,0.028145,0.0696183,0.0703166,0.0298381,-0.0281307,-0.0696326,-0.0703309,-0.0298238,0.028145,0.0696183,0.0703166,0.0298381,-0.0281307,-0.0696326,-0.0703309,-0.0298238,0.028145,0.0696183,0.0703166,0.0298381,-0.0281307,-0.0696326,-0.0703309,-0.0298238,0.028145,0.0696183,0.0703166,0.0298381,-0.0281307,-0.0696326,-0.0703309,-0.0298238,0.028145,0.0696183,0.0703166,0.0298381,-0.0281307,-0.0696326,-0.0703309,-0.0298238,0.028145,0.0696183,0.0703166,0.0298381,-0.0281307,-0.0696326,-0.0703309,-0.0298238,0.028145,0.0696183,0.0703166,0.0298381,-0.0281307,-0.0696326,-0.0703309,-0.0298238,0.028145,0.0696183,0.0703166,0.0298381,-0.0281307,-0.0696326,-0.0703309,-0.0298238,0.028145,0.0696183,0.0703166,0.0298381,-0.0281307,-0.0696326,-0.0703309,-0.0298238,0.028145,0.0696183,0.0703166,0.0298381,-0.0281307,-0.0696326,-0.0703309,-0.0298238,0.028145,0.0696183,0.0703166,0.0298381,-0.0281307,-0.0696326,-0.0703309,-0.0298238,0.028145,0.0696183,0.0703166,0.0298381,-0.0281307,-0.0696326,-0.0703309,-0.0298238,0.028145,0.0696183,0.0703166,0.0298381,-0.0281307,-0.0696326,-0.0703309,-0.0298238,0.028145,0.0696183,0.0703166,0.0298381,-0.0281307,-0.0696326,-0.222914,-0.092816,0.0916389,0.2224,0.222893,0.0928396,-0.0916174,-0.222424,-0.222914,-0.092816,0.0916389,0.2224,0.222893,0.0928396,-0.0916174,-0.222424,-0.222914,-0.092816,0.0916389,0.2224,0.222893,0.0928396,-0.0916174,-0.222424,-0.222914,-0.092816,0.0916389,0.2224,0.222893,0.0928396,-0.0916174,-0.222424,-0.222914,-0.092816,0.0916389,0.2224,0.222893,0.0928396,-0.0916174,-0.222424,-0.222914,-0.092816,0.0916389,0.2224,0.222893,0.0928396,-0.0916174,-0.222424,-0.222914,-0.092816,0.0916389,0.2224,0.222893,0.0928396,-0.0916174,-0.222424,-0.222914,-0.092816,0.0916389,0.2224,0.222893,0.0928396,-0.0916174,-0.222424,-0.222914,-0.092816,0.0916389,0.2224,0.222893,0.0928396,-0.0916174,-0.222424,-0.222914,-0.092816,0.0916389,0.2224,0.222893,0.0928396,-0.0916174,-0.222424,-0.222914,-0.092816,0.0916389,0.2224,0.222893,0.0928396,-0.0916174,-0.222424,-0.222914,-0.092816,0.0916389,0.2224,0.222893,0.0928396,-0.0916174,-0.222424,-0.222914,-0.092816,0.0916389,0.2224,0.222893,0.0928396,-0.0916174,-0.222424,-0.222914,-0.092816,0.0916389,0.2224,0.222893,0.0928396,-0.0916174,-0.222424,-0.222914,-0.092816,0.0916389,0.2224,0.222893,0.0928396,-0.0916174,-0.222424,-0.222914,-0.092816,0.0916389,0.2224,0.222893,0.0928396,-0.0916174,-0.222424,-0.222914,-0.092816,0.0916389,0.2224,0.222893,0.0928396,-0.0916174,-0.222424,-0.222914,-0.092816,0.0916389,0.2224,0.222893,0.0928396,-0.0916174,-0.222424,-0.222914,-0.092816,0.0916389,0.2224,0.222893,0.0928396,-0.0916174,-0.222424,-0.222914,-0.092816,0.0916389,0.2224,0.222893,0.0928396,-0.0916174,-0.222424,-0.222914,-0.092816,0.0916389,0.2224,0.222893,0.0928396,-0.0916174,-0.222424,-0.222914,-0.092816,0.0916389,0.2224,0.222893,0.0928396,-0.0916174,-0.222424,-0.222914,-0.092816,0.0916389,0.2224,0.222893,0.0928396,-0.0916174,-0.222424,-0.222914,-0.092816,0.0916389,0.2224,0.222893,0.0928396,-0.0916174,-0.222424,-0.222914,-0.092816,0.0916389,0.2224,0.222893,0.0928396,-0.0916174,-0.222424,-0.222914,-0.092816,0.0916389,0.2224,0.222893,0.0928396,-0.0916174,-0.222424,-0.222914,-0.092816,0.0916389,0.2224,0.222893,0.0928396,-0.0916174,-0.222424,-0.222914,-0.092816,0.0916389,0.2224,0.222893,0.0928396,-0.0916174,-0.222424,-0.222914,-0.092816,0.0916389,0.2224,0.222893,0.0928396,-0.0916174,-0.222424,-0.222914,-0.092816,0.0916389,0.2224,0.222893,0.0928396,-0.0916174,-0.222424,-0.222914,-0.092816,0.0916389,0.2224,0.222893,0.0928396,-0.0916174,-0.222424,-0.222914,-0.092816,0.0916389,0.2224,0.222893,0.0928396,-0.0916174,-0.222424,-0.0248657,-0.0105443,0.00995076,0.0246138,0.0248607,0.0105494,-0.00994569,-0.0246189,-0.0248657,-0.0105443,0.00995076,0.0246138,0.0248607,0.0105494,-0.00994569,-0.0246189,-0.0248657,-0.0105443,0.00995076,0.0246138,0.0248607,0.0105494,-0.00994569,-0.0246189,-0.0248657,-0.0105443,0.00995076,0.0246138,0.0248607,0.0105494,-0.00994569,-0.0246189,-0.0248657,-0.0105443,0.00995076,0.0246138,0.0248607,0.0105494,-0.00994569,-0.0246189,-0.0248657,-0.0105443,0.00995076,0.0246138,0.0248607,0.0105494,-0.00994569,-0.0246189,-0.0248657,-0.0105443,0.00995076,0.0246138,0.0248607,0.0105494,-0.00994569,-0.0246189,-0.0248657,-0.0105443,0.00995076,0.0246138,0.0248607,0.0105494,-0.00994569,-0.0246189,-0.0248657,-0.0105443,0.00995076,0.0246138,0.0248607,0.0105494,-0.00994569,-0.0246189,-0.0248657,-0.0105443,0.00995076,0.0246138,0.0248607,0.0105494,-0.00994569,-0.0246189,-0.0248657,-0.0105443,0.00995076,0.0246138,0.0248607,0.0105494,-0.00994569,-0.0246189,-0.0248657,-0.0105443,0.00995076,0.0246138,0.0248607,0.0105494,-0.00994569,-0.0246189,-0.0248657,-0.0105443,0.00995076,0.0246138,0.0248607,0.0105494,-0.00994569,-0.0246189,-0.0248657,-0.0105443,0.00995076,0.0246138,0.0248607,0.0105494,-0.00994569,-0.0246189,-0.0248657,-0.0105443,0.00995076,0.0246138,0.0248607,0.0105494,-0.00994569,-0.0246189,-0.0248657,-0.0105443,0.00995076,0.0246138,0.0248607,0.0105494,-0.00994569,-0.0246189,-0.0248657,-0.0105443,0.00995076,0.0246138,0.0248607,0.0105494,-0.00994569,-0.0246189,-0.0248657,-0.0105443,0.00995076,0.0246138,0.0248607,0.0105494,-0.00994569,-0.0246189,-0.0248657,-0.0105443,0.00995076,0.0246138,0.0248607,0.0105494,-0.00994569,-0.0246189,-0.0248657,-0.0105443,0.00995076,0.0246138,0.0248607,0.0105494,-0.00994569,-0.0246189,-0.0248657,-0.0105443,0.00995076,0.0246138,0.0248607,0.0105494,-0.00994569,-0.0246189,-0.0248657,-0.0105443,0.00995076,0.0246138,0.0248607,0.0105494,-0.00994569,-0.0246189,-0.0248657,-0.0105443,0.00995076,0.0246138,0.0248607,0.0105494,-0.00994569,-0.0246189,-0.0248657,-0.0105443,0.00995076,0.0246138,0.0248607,0.0105494,-0.00994569,-0.0246189,-0.0248657,-0.0105443,0.00995076,0.0246138,0.0248607,0.0105494,-0.00994569,-0.0246189,-0.0248657,-0.0105443,0.00995076,0.0246138,0.0248607,0.0105494,-0.00994569,-0.0246189,-0.0248657,-0.0105443,0.00995076,0.0246138,0.0248607,0.0105494,-0.00994569,-0.0246189,-0.0248657,-0.0105443,0.00995076,0.0246138,0.0248607,0.0105494,-0.00994569,-0.0246189,-0.0248657,-0.0105443,0.00995076,0.0246138,0.0248607,0.0105494,-0.00994569,-0.0246189,-0.0248657,-0.0105443,0.00995076,0.0246138,0.0248607,0.0105494,-0.00994569,-0.0246189,-0.0248657,-0.0105443,0.00995076,0.0246138,0.0248607,0.0105494,-0.00994569,-0.0246189,-0.0248657,-0.0105443,0.00995076,0.0246138,0.0248607,0.0105494,-0.00994569,-0.0246189,-0.254572,-0.107442,0.102695,0.252743,0.254687,0.107328,-0.10281,-0.252629,-0.254572,-0.107442,0.102695,0.252743,0.254687,0.107328,-0.10281,-0.252629,-0.254572,-0.107442,0.102695,0.252743,0.254687,0.107328,-0.10281,-0.252629,-0.254572,-0.107442,0.102695,0.252743,0.254687,0.107328,-0.10281,-0.252629,-0.254572,-0.107442,0.102695,0.252743,0.254687,0.107328,-0.10281,-0.252629,-0.254572,-0.107442,0.102695,0.252743,0.254687,0.107328,-0.10281,-0.252629,-0.254572,-0.107442,0.102695,0.252743,0.254687,0.107328,-0.10281,-0.252629,-0.254572,-0.107442,0.102695,0.252743,0.254687,0.107328,-0.10281,-0.252629,-0.254572,-0.107442,0.102695,0.252743,0.254687,0.107328,-0.10281,-0.252629,-0.254572,-0.107442,0.102695,0.252743,0.254687,0.107328,-0.10281,-0.252629,-0.254572,-0.107442,0.102695,0.252743,0.254687,0.107328,-0.10281,-0.252629,-0.254572,-0.107442,0.102695,0.252743,0.254687,0.107328,-0.10281,-0.252629,-0.254572,-0.107442,0.102695,0.252743,0.254687,0.107328,-0.10281,-0.252629,-0.254572,-0.107442,0.102695,0.252743,0.254687,0.107328,-0.10281,-0.252629,-0.254572,-0.107442,0.102695,0.252743,0.254687,0.107328,-0.10281,-0.252629,-0.254572,-0.107442,0.102695,0.252743,0.254687,0.107328,-0.10281,-0.252629,-0.254572,-0.107442,0.102695,0.252743,0.254687,0.107328,-0.10281,-0.252629,-0.254572,-0.107442,0.102695,0.252743,0.254687,0.107328,-0.10281,-0.252629,-0.254572,-0.107442,0.102695,0.252743,0.254687,0.107328,-0.10281,-0.252629,-0.254572,-0.107442,0.102695,0.252743,0.254687,0.107328,-0.10281,-0.252629,-0.254572,-0.107442,0.102695,0.252743,0.254687,0.107328,-0.10281,-0.252629,-0.254572,-0.107442,0.102695,0.252743,0.254687,0.107328,-0.10281,-0.252629,-0.254572,-0.107442,0.102695,0.252743,0.254687,0.107328,-0.10281,-0.252629,-0.254572,-0.107442,0.102695,0.252743,0.254687,0.107328,-0.10281,-0.252629,-0.254572,-0.107442,0.102695,0.252743,0.254687,0.107328,-0.10281,-0.252629,-0.254572,-0.107442,0.102695,0.252743,0.254687,0.107328,-0.10281,-0.252629,-0.254572,-0.107442,0.102695,0.252743,0.254687,0.107328,-0.10281,-0.252629,-0.254572,-0.107442,0.102695,0.252743,0.254687,0.107328,-0.10281,-0.252629,-0.254572,-0.107442,0.102695,0.252743,0.254687,0.107328,-0.10281,-0.252629,-0.254572,-0.107442,0.102695,0.252743,0.254687,0.107328,-0.10281,-0.252629,-0.254572,-0.107442,0.102695,0.252743,0.254687,0.107328,-0.10281,-0.252629,-0.254572,-0.107442,0.102695,0.252743,0.254687,0.107328,-0.10281,-0.252629};
            for (unsigned i = 0; i < _nOut * NUMCOMPS; i++) {
                _rhs[i] = rhs[i];
            }
        }

        virtual void TearDown() {
            free(_Uin);
            free(_rhs);
            free(_rhs_out);
        }

        double _tstop;
        int _nx, _maxstep, _outputInterval;

        unsigned _nOut, _nIn, _nVec, _nAve, _nAve2, _nFlux;

        double* _Uin;
        double* _rhs;
        double* _rhs_out;

        double _retval;
        double _retval_out;
    };

    TEST_F(EulerTest, StepFxn) {
        SetUp({""});

        _retval_out = euler_step(_Uin, _rhs_out);

        ASSERT_LT(Compare(&_retval_out, &_retval, 1), 0);
        ASSERT_LT(Compare(_rhs_out, _rhs, _nOut * NUMCOMPS), 0);

        int stop = 1;
    }
}