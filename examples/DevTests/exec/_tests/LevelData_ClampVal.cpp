#include "Proto.H"
#include "InputParser.H"

using namespace Proto;

PROTO_KERNEL_START
void f_sinF(Point& a_pt, Var<double, 2>& a_data, std::array<double, DIM> a_dx)
{
    std::array<double, DIM> x;
    double v = 1.0;
    for (int dir = 0; dir < DIM; dir++)
    {
        x[dir] = a_pt[dir]*a_dx[dir] + a_dx[dir]/2.0;
        v *= sin(x[dir]);
    }
    a_data(0) = v;
    a_data(1) = v;
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
    int boxSize = 16;
    std::array<bool, DIM> periodicity;
    for (int dir = 0; dir < DIM; dir++) { periodicity[dir] = true; }

    InputArgs args;
    args.add("domainSize",      domainSize);
    args.add("physDomainSize",  physDomainSize);
    args.add("boxSize",         boxSize);
    args.add("periodic_x",      periodicity[0]);
    args.add("periodic_y",      periodicity[1]);
#if DIM > 2
    args.add("periodic_z",      periodicity[2]);
#endif
    args.parse(argc, argv);
    args.print();

    std::array<double, DIM> dx;
    dx.fill(physDomainSize / domainSize);

    Box domainBox = Box::Cube(domainSize);
    Point boxSizeVect = Point::Ones(boxSize);
    ProblemDomain domain(domainBox, periodicity);
    DisjointBoxLayout layout(domain, boxSizeVect);
    LevelBoxData<double, 2> data(layout, Point::Ones());
     
    data.initialize(f_sin, dx);
    h5.writeLevel(dx, data, "DATA_0");
    data.clampVal(0, 0.5, 0);
    h5.writeLevel(dx, data, "DATA_1");
    data.clampVal(-0.5, 0, 1);
    h5.writeLevel(dx, data, "DATA_2");

    LevelBoxData<int> intData(layout, Point::Ones());
    intData.setVal(1);
    h5.writeLevel(intData, "INT_DATA_0");
    intData.clampVal(0,1);
    h5.writeLevel(intData, "INT_DATA_1");


    #ifdef PR_MPI
    MPI_Finalize();
    #endif
    return 0;
}

