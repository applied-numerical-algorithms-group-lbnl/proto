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
    InputArgs args;
    args.parse();
    args.print();

    int domainSize = 64;
    int boxSize = 16;
    int numIter = 3;
    int ghostSize = 1;
    double k = 1;
    double physDomainSize = 1;
    int refRatio = PR_AMR_REFRATIO;
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

