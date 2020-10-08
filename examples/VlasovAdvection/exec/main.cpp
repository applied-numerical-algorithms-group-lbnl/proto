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
#include "VlasovAdvection.H"
#include "Proto_WriteBoxData.H"
#include "Proto_Timer.H"

using namespace std;
using namespace Proto;

int main(int argc, char* argv[])
{
#ifdef PR_MPI
    MPI_Init(&argc,&argv);
#endif

    int pid=procID();

    if(pid==0) {
        cout << "Testing!" << endl;
    }

    LevelBoxData<double,1> rhs;
    LevelBoxData<double,1> phi;
    computeRHS(rhs,phi,1.0);

#ifdef PR_MPI
    MPI_Finalize();
#endif
}
