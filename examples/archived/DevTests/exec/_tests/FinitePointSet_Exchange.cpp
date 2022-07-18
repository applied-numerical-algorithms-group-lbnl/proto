#include "ProtoAMR.H"
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
    std::array<bool, DIM> periodicity;
    for (int dir = 0; dir < DIM; dir++) { periodicity[dir] = true; }

    InputArgs args;
    args.add("domainSize",      domainSize);
    args.add("periodic_x",      periodicity[0]);
    args.add("periodic_y",      periodicity[1]);
#if DIM > 2
    args.add("periodic_z",      periodicity[2]);
#endif
    args.parse(argc, argv);
    args.print();

    Box domainBox = Box::Cube(domainSize);
    ProblemDomain domain(domainBox, periodicity);
    FinitePointSet pointSet(domain);

    if (procID() == 0)
    {
        pointSet.add(Point::Ones(procID()));
    }

    auto points = pointSet.points();
    for (auto iter = points.begin(); iter != points.end(); ++iter)
    {
        pout() << *iter << std::endl;
    }

    #ifdef PR_MPI
    MPI_Finalize();
    #endif
    return 0;
}

