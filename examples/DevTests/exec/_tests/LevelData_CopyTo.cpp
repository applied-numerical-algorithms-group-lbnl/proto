#include "Proto.H"
#include "InputParser.H"

using namespace Proto;

void f_ramp(Point& a_pt, Var<double>& a_data, double a_dx)
{
    a_data(0) = a_pt.sum();
}

int main(int argc, char** argv)
{
    #ifdef PR_MPI
    MPI_Init(&argc, &argv);
    #endif

    // SETUP
    HDF5Handler h5;

    int domainSize = 64;
    int boxSize = 16;
    std::array<bool, DIM> periodicity;
    periodicity.fill(true);

    InputArgs args;
    args.add("domainSize", domainSize);
    args.add("boxSize",    boxSize);
    args.add("periodic_x", periodicity[0]);
    args.add("periodic_y", periodicity[1]);
    args.parse(argc, argv);
    args.print();
    
    double physDomainSize = 1;
    double dx = physDomainSize / domainSize;

    Point boxSizeV = Point::Ones(boxSize);
    Box domainBox = Box::Cube(domainSize);
    ProblemDomain domain(domainBox, periodicity);
    DisjointBoxLayout srcLayout(domain, boxSizeV / 2);
    DisjointBoxLayout dstLayout(domain, boxSizeV);

    LevelBoxData<double> src(srcLayout,  Point::Zeros());
    LevelBoxData<double> sln(dstLayout,  Point::Zeros());
    LevelBoxData<double> dst(dstLayout,  Point::Ones());

    src.initialize(f_ramp, dx);
    sln.initialize(f_ramp, dx);
    dst.setToZero();
    h5.writeLevel(dx, dst, "DST_0");
    src.copyTo(dst);
    h5.writeLevel(dx, src, "SRC");
    h5.writeLevel(dx, dst, "DST_1");

    dst.increment(sln, -1);
    std::cout << "Error: " << dst.absMax() << std::endl;

    #ifdef PR_MPI
    MPI_Finalize();
    #endif
    return 0;
}

