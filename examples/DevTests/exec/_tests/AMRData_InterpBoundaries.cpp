#include "Proto.H"
#include "InputParser.H"
#include "func.H"

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
    int numIter = 3;
    int numLevels = 2;
    int ghostSize = 1;
    int order = 4;
    int refRatio = 4;
    std::array<bool, DIM> periodicity;
    for (int dir = 0; dir < DIM; dir++) { periodicity[dir] = true; }

    InputArgs args;
    args.add("domainSize", domainSize);
    args.add("boxSize",    boxSize);
    args.add("ghostSize",  ghostSize);
    args.add("numIter",    numIter);
    args.add("numLevels",  numLevels);
    args.add("periodic_x", periodicity[0]);
    args.add("periodic_y", periodicity[1]);
    args.add("order",      order);
    args.parse(argc, argv);
    args.print();
    
    double k = 1;
    double physDomainSize = 1;

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

        AMRData<double> Phi(grid,            Point::Ones(ghostSize));
        AMRData<double> PhiSln(grid,         Point::Ones(ghostSize));
        LevelBoxData<double> PhiErr(grid[1], Point::Ones(ghostSize));
        
        Phi.initConvolve(dx, f_wave, k);
        PhiSln.initConvolve(dx, f_wave, k);
        for (auto iter = layouts[1].begin(); iter.ok(); ++iter)
        {
            auto& phi_i = Phi[1][*iter];
            auto& sln_i = PhiSln[1][*iter];
            forallInPlace_p(f_clearExterior, phi_i, sln_i, refinedRegion);
        }
        PhiErr.setToZero();
        h5.writeAMRData(dx, Phi, "Phi_N%i_0", nn);
        h5.writeAMRData(dx, PhiSln, "PhiSln_N%i", nn);
        
        auto INTERP = InterpStencil<double>::Build(order, refRatio);
        //auto INTERP = InterpStencil<double>::Quadratic(refRatio);
        interpBoundaries(Phi[0], Phi[1], INTERP);
        h5.writeAMRData(dx, Phi, "Phi_N%i_1", nn);

        Reduction<double, Abs> errNorm;
        for (auto iter = layouts[1].begin(); iter.ok(); ++iter)
        {
            auto& phi_i = Phi[1][*iter];
            auto& sln_i = PhiSln[1][*iter];
            auto& err_i = PhiErr[*iter];

            phi_i.copyTo(err_i);
            err_i -= sln_i;
            forallInPlace_p(f_clearInterior, err_i, err_i, refinedRegion);
            err_i.reduce(errNorm);
        }
        h5.writeLevel(dx/PR_AMR_REFRATIO, PhiErr, "PhiErr_N%i", nn);
        err[nn] = errNorm.fetch();
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

