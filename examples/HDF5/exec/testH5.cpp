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
    
    a_data(0) = sin(2.0*M_PI*(x[0] + x[1]));
    a_data(1) = sin(2.0*M_PI*(x[0] - x[1]));
    a_data(2) = a_data(0) + a_data(1);
}

int main(int argc, char** argv)
{
    double L = 1.0;
    int domainSize = 16;
    Point boxSize = Point::Ones(8);
    Box domainBox = Box::Cube(domainSize);
    array<bool, DIM> periodicity;
    for (int ii = 0; ii < DIM; ii++){ periodicity[ii] = true; }
    ProblemDomain domain(domainBox, periodicity);
    DisjointBoxLayout layout(domain, boxSize);
    
    LevelBoxData<double, D> writtenData(layout, Proto::Point::Ones(2)); 
    LevelBoxData<double, D> readData(layout, Proto::Point::Ones(2)); 
    array<double, DIM> dx;
    for (int ii = 0; ii < DIM; ii++){ dx[ii] = L/domainSize; }

    for (auto iter = writtenData.begin(); *iter != iter.end(); ++iter)
    {
        auto& patch = writtenData[*iter];
        forallInPlace_p(f_init, patch, dx);
    }

    HDF5Handler h5;
    h5.writeLevel(writtenData, {"alpha", "beta", "gamma"}, dx,"LevelData");
    h5.readLevel(readData, "LevelData");

    double error = 0.0;
    for (auto iter = writtenData.begin(); *iter != iter.end(); ++iter)
    {
        auto& wpatch = writtenData[*iter];
        auto& rpatch = readData[*iter];
        wpatch -= rpatch;
        error = max(error, wpatch.absMax());
    }
    std::cout << "Write -> Read Error: " << error << std::endl;
}
