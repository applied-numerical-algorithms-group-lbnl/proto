#include "Proto.H"
#include "InputParser.H"
#include "func.H"

using namespace Proto;

void f_foo(Point& a_pt, Var<double> a_data, double a_dx, Box a_box)
{
    if (a_box.contains(a_pt))
    {
        a_data(0) = -1;
    } else {
        a_data(0) = 1;
    }
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
    int numLevels = 2;
    int ghostSize = 1;
    int order = 4;
    double k = 1;
    double physDomainSize = 1;
    int refRatio = PR_AMR_REFRATIO;
    std::array<bool, DIM> periodicity;
    for (int dir = 0; dir < DIM; dir++) { periodicity[dir] = true; }

    args.set("domainSize", &domainSize);
    args.set("boxSize",    &boxSize);
    args.set("ghostSize",  &ghostSize);
    args.set("numIter",    &numIter);
    args.set("numLevels",  &numLevels);
    args.set("periodic_x", &periodicity[0]);
    args.set("periodic_y", &periodicity[1]);
    args.set("order",      &order);

    double err[numIter];
    for (int nn = 0; nn < numIter; nn++)
    {
        // BUILD GRIDS
        double dx = physDomainSize / domainSize;

        std::vector<DisjointBoxLayout> layouts;
        layouts.resize(numLevels);
                
        // coarse layout
        Point boxSizeV = Point::Ones(boxSize);
        Box domainBox = Box::Cube(domainSize);
        ProblemDomain domain(domainBox, periodicity);
        layouts[0].define(domain, boxSizeV);

        Box refinedRegion = domainBox;
        for (int lvl = 1; lvl < numLevels; lvl++)
        {
            Point prevSize = refinedRegion.high() - refinedRegion.low() + Point::Ones();
            refinedRegion = refinedRegion.grow(-prevSize / 4).refine(refRatio);
            
            Box refinedRegionPatches = refinedRegion.coarsen(boxSizeV);
            std::vector<Point> fineLayoutPatches;
            bool skip = true;
            for (auto iter = refinedRegionPatches.begin(); iter.ok(); ++iter)
            {
                if (skip)
                {
                    skip = false;
                    continue;
                }
                fineLayoutPatches.push_back(*iter);
            }
            ProblemDomain fineDomain = layouts[lvl-1].domain().refine(Point::Ones(refRatio));
            layouts[lvl].define(fineDomain, fineLayoutPatches, boxSizeV);
        }

        AMRGrid grid(layouts, numLevels);

        AMRData<double> Phi(grid,            Point::Zeros());
        
        Phi.initialize(dx, f_foo, refinedRegion);
        h5.writeAMRData(dx, Phi, "Phi_N%i_0", nn);
        double integral = Phi.integrateAbs(dx);
        std::cout << "Integral is " << integral << std::endl;
        
        err[nn] = (integral - 1);
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

