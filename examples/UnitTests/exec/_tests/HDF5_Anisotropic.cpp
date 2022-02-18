#include "Proto.H"
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
    InputArgs args;
    args.parse();
    args.print();

    int domainSize = 64;
    int boxSize = 16;
    int numIter = 3;
    int ghostSize = 1;
    double k = 1;
    double physDomainSize = 1;
    int refRatio = 4;
    int numLevels = 2;
    std::array<bool, DIM> periodicity;
    for (int dir = 0; dir < DIM; dir++) { periodicity[dir] = true; }

    args.set("domainSize", &domainSize);
    args.set("boxSize",    &boxSize);
    args.set("numLevels",  &numLevels);
    args.set("numIter",    &numIter);
    args.set("periodic_x", &periodicity[0]);
    args.set("periodic_y", &periodicity[1]);

    double err[numIter];
    for (int nn = 0; nn < numIter; nn++)
    {
        // BUILD GRIDS
        double dx = physDomainSize / domainSize;
        std::array<double, DIM> dxVect;
        dxVect[0] = dx;
        dxVect[1] = dx;

        std::vector<Point> refRatios;
        Point refRatioVect(1,2,2,2,2,2);
        refRatios.resize(numLevels-1, refRatioVect);
        //refRatios.resize(numLevels-1, Point::Ones(refRatio));
        std::vector<DisjointBoxLayout> layouts;
        layouts.resize(numLevels);
                
        Point boxSizeV = Point::Ones(boxSize);
        Box domainBox = Box::Cube(domainSize);
        ProblemDomain domain(domainBox, periodicity);
        layouts[0].define(domain, boxSizeV);

        Box refinedRegion = domainBox;
        for (int lvl = 1; lvl < numLevels; lvl++)
        {
            Point prevSize = refinedRegion.high() - refinedRegion.low() + Point::Ones();
            refinedRegion = refinedRegion.grow(-prevSize / 4).refine(refRatios[lvl-1]);
            Box refinedRegionPatches = refinedRegion.coarsen(boxSizeV);
            std::vector<Point> fineLayoutPatches;
            for (auto iter = refinedRegionPatches.begin(); iter.ok(); ++iter)
            {
                fineLayoutPatches.push_back(*iter);
            }
            ProblemDomain fineDomain = layouts[lvl-1].domain().refine(refRatios[lvl-1]);
            layouts[lvl].define(fineDomain, fineLayoutPatches, boxSizeV);
        }
        AMRGrid grid(layouts, refRatios, numLevels);
        for (int lvl = 0; lvl < numLevels; lvl++)
        {
            grid[lvl].print();
        }
        AMRData<double> data(grid, Point::Ones());
        data.setToZero();
        h5.writeAMRData(dxVect, data, "DATA");

        

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

