#include "Proto.H"
#include "Proto_LevelBoundaryRegister.H"
#include "Proto_LevelBC.H"
#include "InputParser.H"

using namespace Proto;

unsigned int foo(Point& a_pt, Var<double, 1>& a_data, int a_dir)
{
    a_data(0) = a_pt[a_dir];
    return 0;
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
    args.set("periodic_x", &periodicity[0]);
    args.set("periodic_y", &periodicity[1]);

    Point boxSizeV = Point::Ones(boxSize);
    Box domainBox = Box::Cube(domainSize);
    ProblemDomain crseDomain(domainBox, periodicity);
    DisjointBoxLayout crseLayout(crseDomain, boxSizeV);

    LevelBC<NEUMANN, double, 1> bc0(crseLayout, Point::Ones());
    bc0.set(0, Side::Lo, 0, 10);
    bc0.set(0, Side::Hi, 0, 11);
    bc0.set(1, Side::Lo, 0, 20);
    bc0.set(1, Side::Hi, 0, 21);
   
    int x = 0; int y = 1;
    LevelBC<NEUMANN, double, 1> bc1(crseLayout, Point::Ones());
    bc1.set(0, Side::Lo, foo, x);
    bc1.set(0, Side::Hi, foo, x);
    bc1.set(1, Side::Lo, foo, y);
    bc1.set(1, Side::Hi, foo, y);
    
    LevelBoxData<double, 1> data(crseLayout, Point::Ones());
    data.setToZero();

    h5.writeLevel(data, "DATA_0");
    for (auto iter = data.begin(); iter.ok(); ++iter)
    {
        auto& data_i = data[*iter];
        bc0.apply(data_i, *iter);
    }
    h5.writeLevel(data, "DATA_1");
    
    data.setToZero();
    h5.writeLevel(data, "DATA_2");
    for (auto iter = data.begin(); iter.ok(); ++iter)
    {
        auto& data_i = data[*iter];
        bc1.apply(data_i, *iter);
    }
    h5.writeLevel(data, "DATA_3");

    #ifdef PR_MPI
    MPI_Finalize();
    #endif
    return 0;
}

