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
    int boxSize = 16;
    std::array<bool, DIM> periodicity;
    for (int dir = 0; dir < DIM; dir++) { periodicity[dir] = true; }
    int testNum = 0;

    InputArgs args;
    args.add("testNum",    testNum);
    args.add("domainSize", domainSize);
    args.add("boxSize",    boxSize);
    args.add("periodic_x", periodicity[0]);
    args.add("periodic_y", periodicity[1]);
    
    args.parse(argc, argv);
    args.print();

    if (testNum == 0)
    {
        pout() << "Testing DisjointBoxLayout::onDomainBoundary" << std::endl;
        Point boxSizeVect = Point::Ones(boxSize);
        Box domainBox = Box::Cube(domainSize);
        ProblemDomain domain(domainBox, periodicity);
        DisjointBoxLayout layout(domain, boxSizeVect);

        layout.print();
        std::vector<Point> points = {Point::Zeros(), Point::Basis(0), Point::Basis(1), Point::Ones()};
        for (auto iter = points.begin(); iter != points.end(); ++iter)
        {
            pout() << "onDomainBoundary(" << *iter << "): " << layout.onDomainBoundary(*iter) << std::endl;
        }
    } else if (testNum == 1)
    {
        pout() << "Testing DisjointBoxLayout::onDomainBoundary" << std::endl;
        Point boxSizeVect = Point::Ones(boxSize);
        Box domainBox = Box::Cube(domainSize);
        Box patches = (domainBox.shift(Point::Basis(0, boxSize)) & domainBox).coarsen(boxSizeVect);
        std::vector<Point> patchPoints;
        for (auto biter = patches.begin(); biter != patches.end(); ++biter)
        {
            patchPoints.push_back(*biter);
        }
        ProblemDomain domain(domainBox, periodicity);
        DisjointBoxLayout layout(domain, patchPoints, boxSizeVect);

        layout.print();
        for (auto iter = patchPoints.begin(); iter != patchPoints.end(); ++iter)
        {
            pout() << "onLevelBoundary(" << *iter << "): " << layout.onLevelBoundary(*iter) << std::endl;
        }

    }



    #ifdef PR_MPI
    MPI_Finalize();
    #endif
    return 0;
}

