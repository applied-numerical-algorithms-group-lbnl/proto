#include "Proto.H"
#include "InputParser.H"

using namespace Proto;

PROTO_KERNEL_START
void f_ramp_0(Point& a_pt, Var<double>& a_data, double a_dx)
{
    a_data(0) = a_pt.sum();
}
PROTO_KERNEL_END(f_ramp_0, f_ramp);

int main(int argc, char** argv)
{
    #ifdef PR_MPI
    MPI_Init(&argc, &argv);
    #endif

    // SETUP
    HDF5Handler h5;

    int domainSize = 64;
    int boxSize = 16;
    int testNum = 0;
    std::array<bool, DIM> periodicity;
    periodicity.fill(true);

    InputArgs args;
    args.add("testNum",    testNum);
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

    pout() << "POINT 0" << std::endl;    
    DisjointBoxLayout srcLayout;
    DisjointBoxLayout dstLayout;
    if (testNum == 0)
    {
        srcLayout.define(domain, boxSizeV / 2);
        dstLayout.define(domain, boxSizeV);
    } else if (testNum == 1)
    {
        srcLayout.define(domain, boxSizeV);
        dstLayout.define(domain, boxSizeV / 2);
    } else {
        pout() << "Invalid testNum: " << testNum << std::endl;
        std::abort();
    }
    pout() << "POINT 1" << std::endl;    

    LevelBoxData<double> src(srcLayout,  Point::Zeros());
    LevelBoxData<double> sln(dstLayout,  Point::Zeros());
    LevelBoxData<double> dst(dstLayout,  Point::Ones());

    pout() << "POINT 2" << std::endl;    
    src.initialize(f_ramp, dx);
    sln.initialize(f_ramp, dx);
    dst.setToZero();
    pout() << "POINT 3" << std::endl;    
    //h5.writeLevel(dx, dst, "DST_0");
    pout() << "POINT 4" << std::endl;    
    src.copyTo(dst);
    pout() << "POINT 5" << std::endl;    
    //h5.writeLevel(dx, src, "SRC");
    //h5.writeLevel(dx, dst, "DST_1");

    dst.increment(sln, -1);
    std::cout << "Error: " << dst.absMax() << std::endl;

    #ifdef PR_MPI
    MPI_Finalize();
    #endif
    return 0;
}

