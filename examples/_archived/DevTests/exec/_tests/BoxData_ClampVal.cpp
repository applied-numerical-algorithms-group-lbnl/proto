#include "Proto.H"
#include "InputParser.H"

using namespace Proto;

PROTO_KERNEL_START
void f_sinF(Point& a_pt, Var<double> a_data, std::array<double, DIM> a_dx)
{
    std::array<double, DIM> x;
    a_data(0) = 1.0;
    for (int dir = 0; dir < DIM; dir++)
    {
        x[dir] = a_pt[dir]*a_dx[dir] + a_dx[dir]/2.0;
        a_data(0) *= sin(x[dir]);
    }
}
PROTO_KERNEL_END(f_sinF, f_sin);

int main(int argc, char** argv)
{
    #ifdef PR_MPI
    MPI_Init(&argc, &argv);
    #endif

    // SETUP
    HDF5Handler h5;

    int domainSize = 64;
    double physDomainSize = 2*M_PI;
    double max = 0.5;
    double min = -0.75;
    InputArgs args;
    args.add("domainSize",      domainSize);
    args.add("physDomainSize",  physDomainSize);
    args.add("max", max);
    args.add("min", min);
    args.parse(argc, argv);
    args.print();

    Box domainBox = Box::Cube(domainSize);
    std::array<double, DIM> dx;
    dx.fill(physDomainSize / domainSize);

    auto data = forall_p<double>(f_sin, domainBox, dx);
    
    data.clampVal(min, max);
   
    bool testPassed = true;
    testPassed &= data.max() <= max;
    testPassed &= data.min() >= min;
    std::cout << "Test Passed: " << testPassed << std::endl;

    #ifdef PR_MPI
    MPI_Finalize();
    #endif
    return 0;
}

