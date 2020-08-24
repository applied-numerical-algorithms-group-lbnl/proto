#include <cstdio>
#include <cstring>
#include <cassert>
#include <cmath>

#include <vector>
#include <memory>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "Proto.H"
#include "EulerLevelDataRK4.H"
#include "Proto_WriteBoxData.H"
#include "Proto_Timer.H"

#define PI 3.141592653589793

using namespace std;
using namespace Proto;

int main(int argc, char* argv[])
{
#ifdef PR_MPI
    MPI_Init(&argc,&argv);
#endif

    RK4<EulerLevelDataState,EulerLevelDataRK4Op,EulerLevelDataDX> rk4;

#ifdef PR_MPI
    MPI_Finalize();
#endif
}
