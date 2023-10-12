#include <iostream>
#include "Proto.H"

#define D 3

using namespace Proto;

void f_init(Point& a_pt, Var<double, D>& a_data, array<double, DIM> a_dx)
{
   
    double x[DIM];
   
    for (int ii = 0; ii < DIM; ii++)
    {
        x[ii] = a_pt[ii]*a_dx[ii] + a_dx[ii]/2.0;
    }
    
    #if DIM == 2
    a_data(0) = sin(2.0*M_PI*(x[0] + x[1]));
    a_data(1) = sin(2.0*M_PI*(x[0] - x[1]));
    #elif DIM == 3
    a_data(0) = sin(2.0*M_PI*(x[0] + x[1] + x[2]));
    a_data(1) = sin(2.0*M_PI*(x[0] - x[1] + x[2]));
    #endif
    a_data(2) = a_data(0) + a_data(1);
}

int main(int argc, char** argv)
{
    PR_TIMER_SETFILE("testH5.timer");
    PR_TIME("main");
    #ifdef PR_MPI
    MPI_Init(&argc, &argv);
    #endif
    
    double L = 1.0;
    int domainSize = 16;
    if (argc > 1 )
    {
        domainSize = atoi(argv[1]);    
    }

    Point boxSize = Point::Ones(domainSize / 4);
    Box domainBox = Box::Cube(domainSize);
    array<bool, DIM> periodicity;
    for (int ii = 0; ii < DIM; ii++){ periodicity[ii] = false; }
    ProblemDomain domain(domainBox, periodicity);
    DisjointBoxLayout layout(domain, boxSize);
    
    LevelBoxData<double, D> writtenData(layout, Proto::Point::Ones(2)); 
    LevelBoxData<double, D> readData(layout, Proto::Point::Ones(2)); 
    array<double, DIM> dx;
    for (int ii = 0; ii < DIM; ii++){ dx[ii] = L/domainSize; }

    for (auto iter = writtenData.begin(); iter != iter.end(); ++iter)
    {
        auto& patch = writtenData[*iter];
        forallInPlace_p(f_init, patch, dx);
    }

    HDF5Handler h5;
    /*
    for (auto iter = writtenData.begin(); iter != iter.end() ;++iter)
    {
        auto& patch = writtenData[*iter];
        h5.writePatch(patch, {"alpha", "beta", "gamma"}, dx, "PatchData_%i", (*iter).intIndex());
    }
    */
    h5.writeLevel({"alpha", "beta", "gamma"}, L/domainSize, writtenData, "LevelData");
    barrier();
    h5.readLevel(readData, "LevelData");
    double error = 0.0;
    for (auto iter = writtenData.begin(); iter != iter.end(); ++iter)
    {
        auto& wpatch = writtenData[*iter];
        auto& rpatch = readData[*iter];
        wpatch -= rpatch;
        error = max(error, wpatch.absMax());
    }
    std::cout << "Write -> Read Error: " << error << std::endl;
    #ifdef PR_MPI
    MPI_Finalize();
    #endif
    PR_TIMER_REPORT();
}
