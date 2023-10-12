#include "ProtoAMR.H"
#include "InputParser.H"

using namespace Proto;

int main(int argc, char** argv)
{
    #ifdef PR_MPI
    MPI_Init(&argc, &argv);
    #endif

    // SETUP
    using Proto::pout;
    HDF5Handler h5;

    int domainSize = 32;
    int boxSize = 32;
    int numIter = 3;
    int ghostSize = 1;
    double k = 1;
    double physDomainSize = 1;
    int refRatio = 4;
    int numLevels = 2;
    std::array<bool, DIM> periodicity;
    for (int dir = 0; dir < DIM; dir++) { periodicity[dir] = true; }

    InputArgs args;
    args.add("domainSize", domainSize);
    args.add("boxSize",    boxSize);
    args.add("numLevels",  numLevels);
    args.add("numIter",    numIter);
    args.add("periodic_x", periodicity[0]);
    args.add("periodic_y", periodicity[1]);
    args.parse(argc, argv);
    args.print();

    double dx = physDomainSize / domainSize;

    Box domainBox = Box::Cube(domainSize);
    Point boxSizeVect = Point::Ones(boxSize);
    ProblemDomain domain(domainBox, periodicity);
    DisjointBoxLayout layout(domain, boxSizeVect);

    Box refinedRegion = Box::Cube(domainSize / 2).refine(refRatio);
    DisjointBoxLayout fineLayout(domain.refine(refRatio), refinedRegion, boxSizeVect);

    std::vector<DisjointBoxLayout> layouts;
    layouts.push_back(layout);
    layouts.push_back(fineLayout);
    
    std::vector<Point> refRatios;
    refRatios.push_back(Point::Ones(refRatio));

    AMRGrid grid(layouts, refRatios, 2);
    AMRData<double> data(grid, Point::Zeros());
    data.setToZero();
    h5.writeAMRData(dx, data, "DATA");
    
#ifdef PR_MPI
    MPI_Finalize();
#endif
    return 0;
}

