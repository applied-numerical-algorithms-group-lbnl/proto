#include <cstring>
#include <string>
using std::string;
using std::to_string;

#include <gtest/gtest.h>
using namespace testing;

#include <BenchmarkTest.hpp>
#include <util/Lists.hpp>

#define PI 3.141592653589793
#define NGHOST 4

#include "Proto.H"
#include "EulerOp.H"

// Import generated code
//#include "euler_step.h"
#if DIM>2
//#include "euler_step_3d.h"
//#include "euler_step_3d_fuse.h"
//#include "euler_step_3d_unroll.h"
#include "euler_step_3d_opt.h"
#define DATA_FILE "data/Uin_3d.csv"
#else
#include "euler_step_2d.h"
#define DATA_FILE "data/Uin_2d.csv"
#endif

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

            _nOut = 1;
            _nIn = 1;
            _nVec = 1;
//            _nAve = (NUMCELLS + NGHOST/2) * (NUMCELLS + NGHOST + NGHOST/2);
//            _nAve2 = (NUMCELLS + 1) * (NUMCELLS + NGHOST + NGHOST/2);
//            _nFlux = (NUMCELLS + 1) * (NUMCELLS + NGHOST);

            for (unsigned d = 0; d < DIM; d++) {
                _nOut *= NUMCELLS;
                _nIn *= (NUMCELLS + 2 * NGHOST);
                _nVec *= (NUMCELLS + NGHOST);
            }

            _Uin = (double*) calloc(_nIn * NUMCOMPS, sizeof(double));
            vector<double> Uinit(_nIn * NUMCOMPS);
//            fprintf(stderr, "nIn=%u,nOut=%u,nGhost=%u,nCell=%u,nComp=%u,exp=%u,cap=%lu,outsize=%u,file='%s'\n", _nIn,
//                    _nOut, NGHOST, NUMCELLS, NUMCOMPS, _nIn * NUMCOMPS, Uinit.capacity(), _nOut * NUMCOMPS, DATA_FILE);
            Lists::read<double>(Uinit, DATA_FILE);
            std::copy(Uinit.begin(), Uinit.end(), _Uin);

            _velmax = 0.0; // 2D16=2.76852; 3D64=4.0009
            _rhs = (double*) calloc(_nOut * NUMCOMPS, sizeof(double));
            _rhs_out = (double*) calloc(_nOut * NUMCOMPS, sizeof(double));

            // Setup evaluation
#if DIM>2
            _dbx0 = Box(Point(0,0,0), Point(NUMCELLS-1,NUMCELLS-1,NUMCELLS-1));
#else
            _dbx0 = Box(Point(0,0), Point(NUMCELLS-1,NUMCELLS-1));
#endif
            Box dbx = _dbx0.grow(NGHOST);
            _Uave = BoxData<double,NUMCOMPS>(_Uin, dbx);
            _dxdu = BoxData<double,NUMCOMPS>(_dbx0);
        }

        virtual void Execute() {
            _velmax_out = euler_step(_Uin, _rhs_out);
        }

        // Execute reference code for verification
        virtual void Evaluate() {
            _velmax = EulerOp::step(_dxdu, _Uave, _dbx0);
        }

        virtual void Assert() {
            // Copy dx_du data into _rhs:
            memcpy(_rhs, _dxdu.data(), _dxdu.size() * sizeof(double));

            ASSERT_FALSE(isnan(_rhs_out[0]));
            ASSERT_FALSE(isnan(_rhs_out[_nOut * NUMCOMPS - 1]));
            ASSERT_LT(Compare(_rhs_out, _rhs, _nOut * NUMCOMPS), 0);
            ASSERT_FALSE(isnan(_velmax_out));
            ASSERT_LT(Compare(&_velmax_out, &_velmax, 1), 0);
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

        double _velmax;
        double _velmax_out;

        Box _dbx0;
        BoxData<double,NUMCOMPS> _Uave;
        BoxData<double,NUMCOMPS> _dxdu;
    };

    TEST_F(EulerTest, StepFxn) {
        SetUp({""});
        NumRuns(3);
        Run();
        Verify();
        Assert();
    }
}