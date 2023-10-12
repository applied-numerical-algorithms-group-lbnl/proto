#include "Proto.H"
#include "InputParser.H"

using namespace Proto;

void f_ramp(Point& a_pt, Var<int>& a_data)
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

    int domainSize = 32;
    double physDomainSize = 1.0;
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

    double dx = physDomainSize / domainSize;

    Point boxSizeV = Point::Ones(boxSize);
    Box domainBox = Box::Cube(domainSize);
    ProblemDomain domain(domainBox, periodicity);
    DisjointBoxLayout layout(domain, boxSizeV);

    LevelBoxData<int> data(layout, Point::Zeros());
    for (auto iter = layout.begin(); iter.ok(); ++iter)
    {
        forallInPlace_p(f_ramp, data[*iter]);
    }
    
    h5.writeLevel(dx, data, "LEVEL_DATA");

    for (auto iter = layout.begin(); iter.ok(); ++iter)
    {
        h5.writePatch(dx, data[*iter], "PATCH_%i", (*iter).global());
    }

    #ifdef PR_MPI
    MPI_Finalize();
    #endif
    return 0;
}
