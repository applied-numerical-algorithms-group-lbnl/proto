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
    InputArgs args;
    args.parse();
    args.print();

    int domainSize = 64;
    int boxSize = 16;
    int numIter = 3;
    int ghostSize = 1;
    double k = 1;
    double physDomainSize = 1;
    std::array<bool, DIM> periodicity;
    for (int dir = 0; dir < DIM; dir++) { periodicity[dir] = true; }

    args.set("domainSize", &domainSize);
    args.set("boxSize",    &boxSize);
    args.set("numIter",    &numIter);
    args.set("periodic_x", &periodicity[0]);
    args.set("periodic_y", &periodicity[1]);

    double err[numIter];
    for (int nn = 0; nn < numIter; nn++)
    {
        // BUILD GRIDS
        double dx = physDomainSize / domainSize;

        Point boxSizeV = Point::Ones(boxSize);
        Box domainBox = Box::Cube(domainSize);
        ProblemDomain domain(domainBox, periodicity);
        DisjointBoxLayout srcLayout(domain, boxSizeV);

        DisjointBoxLayout dstLayout(domain, boxSizeV / 2);

        LevelBoxData<double> src(dstLayout,  Point::Zeros());
        LevelBoxData<double> dst(srcLayout,  Point::Ones());
        
        src.initialize(f_ramp, dx);
        dst.setToZero();
        h5.writeLevel(dx, dst, "DST_0");
        src.copyTo(dst);
        h5.writeLevel(dx, src, "SRC");
        h5.writeLevel(dx, dst, "DST_1");

        pout() << "Error Rate: " << err[nn] << std::endl;
        domainSize *= 2;
    }
        
    for (int ii = 1; ii < numIter; ii++)
    {
        pout() << "Convergence Rate: " << log(err[ii-1] / err[ii]) / log(2.0) << std::endl;
    }

    #ifdef PR_MPI
    MPI_Finalize();
    #endif
    return 0;
}

