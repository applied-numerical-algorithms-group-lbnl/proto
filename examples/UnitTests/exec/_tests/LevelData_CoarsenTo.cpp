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
        double cdx = dx;
        double fdx = dx / refRatio;

        // coarse layout
        Point boxSizeV = Point::Ones(boxSize);
        Box domainBox = Box::Cube(domainSize);
        ProblemDomain crseDomain(domainBox, periodicity);
        DisjointBoxLayout crseLayout(crseDomain, boxSizeV);

        ProblemDomain fineDomain = crseDomain.refine(Point::Ones(refRatio));
        Box refinedRegion = Box::Cube(boxSize).shift(Point::Ones(boxSize)).refine(Point::Ones(refRatio));
        DisjointBoxLayout fineLayout(fineDomain, refinedRegion, boxSizeV);

        LevelBoxData<double> Fine(fineLayout,  Point::Ones());
        LevelBoxData<double> Crse(crseLayout,  Point::Ones());
        LevelBoxData<double> CF(crseLayout,    Point::Ones(2));
        LevelBoxData<double> Error(crseLayout, Point::Zeros());
        
        Fine.initConvolve(f_wave, fdx, k);
        Crse.initConvolve(f_wave, cdx, k);
        CF.setToZero();
        Fine.coarsenTo(CF, Point::Ones(refRatio));
        
        for (auto iter = crseLayout.begin(); iter.ok(); ++iter)
        {
            CF[*iter].printData();
        }
        
        h5.writeLevel(fdx, Fine, "SRC_N%i", nn);
        h5.writeLevel(cdx, CF, "DST_N%i", nn);

        for (auto iter = crseLayout.begin(); iter.ok(); ++iter)
        {
            auto& error_i = Error[*iter];
            auto& cf_i    = CF[*iter];
            auto& crse_i  = Crse[*iter];

            cf_i.copyTo(error_i);
            error_i -= crse_i;
        }
        
        err[nn] = Error.absMax();

        pout() << "Error Rate: " << err[nn] << std::endl;
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

