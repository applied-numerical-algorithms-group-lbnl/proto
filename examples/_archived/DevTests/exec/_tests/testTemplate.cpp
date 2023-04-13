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
    double physDomainSize = 1.0;
    int boxSize = 16;
    int numIter = 3;
    int refRatio = 2;
    std::array<bool, DIM> periodicity;
    for (int dir = 0; dir < DIM; dir++) { periodicity[dir] = true; }

    InputArgs args;
    args.add("domainSize",      domainSize);
    args.add("physDomainSize",  physDomainSize);
    args.add("boxSize",         boxSize);
    args.add("numIter",         numIter);
    args.add("refRatio",        refRatio);
    args.add("periodic_x",      periodicity[0]);
    args.add("periodic_y",      periodicity[1]);
#if DIM > 2
    args.add("periodic_z",      periodicity[2]);
#endif
    args.parse(argc, argv);
    args.print();

    double err[numIter];
    for (int nn = 0; nn < numIter; nn++)
    {
        // BUILD GRIDS
        double dx = physDomainSize / domainSize;
        
        // coarse layout
        Point boxSizeV = Point::Ones(boxSize);
        Box domainBox = Box::Cube(domainSize);
        ProblemDomain crseDomain(domainBox, periodicity);
        DisjointBoxLayout crseLayout(crseDomain, boxSizeV);

        // fine layout
        Box refinedRegion = Box::Cube(domainSize / 2).shift(Point::Ones(domainSize / 4));
        refinedRegion = refinedRegion.refine(refRatio);
        Box refinedRegionPatches = refinedRegion.coarsen(boxSizeV);
        std::vector<Point> fineLayoutPatches;
        for (auto iter = refinedRegionPatches.begin(); iter.ok(); ++iter)
        {
            fineLayoutPatches.push_back(*iter);
        }
        ProblemDomain fineDomain = crseDomain.refine(Point::Ones(refRatio));
        DisjointBoxLayout fineLayout(fineDomain, fineLayoutPatches, boxSizeV);

        // test code here

        pout() << "Error: " << err[nn] << std::endl;
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

