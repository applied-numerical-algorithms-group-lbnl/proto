#include "Proto.H"
#include "InputParser.H"

using namespace Proto;

int main(int argc, char** argv)
{
    #ifdef PR_MPI
    MPI_Init(&argc, &argv);
    #endif

    // SETUP
    HDF5Handler h5;

    int domainSize = 64;
    double physDomainSize = 2*M_PI;
    InputArgs args;
    args.add("domainSize",      domainSize);
    args.add("physDomainSize",  physDomainSize);
    args.parse(argc, argv);
    args.print();

    Box domainBox = Box::Cube(domainSize);
    BoxData<double, 2> data(domainBox);
    bool testPassed = true;
    data.setToZero();
    testPassed &= (data.max() < 1e-12);
    testPassed &= (data.min() > -1e-12);
    data.setVal(7);
    testPassed &= ((data.max() - 7) < 1e-12);
    testPassed &= ((data.min() - 7) > -1e-12);

    #ifdef PR_MPI
    MPI_Finalize();
    #endif
    return 0;
}

