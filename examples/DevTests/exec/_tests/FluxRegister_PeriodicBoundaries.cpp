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
    int numIter = 1;
    double physDomainSize = 1;
    int refRatio = 2;
    std::array<bool, DIM> periodicity;
    for (int dir = 0; dir < DIM; dir++) { periodicity[dir] = true; }

    args.set("domainSize", &domainSize);
    args.set("boxSize",    &boxSize);
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

        ProblemDomain fineDomain = crseDomain.refine(refRatio);
        DisjointBoxLayout fineLayout(fineDomain, boxSizeV);
        // test code here

        crseLayout.print();
        fineLayout.print();

        LevelBoxData<double> data(crseLayout, Point::Zeros());
        data.setToZero();
        LevelFluxRegister<double> FR(crseLayout, fineLayout, Point::Ones(refRatio));
        FR.print();
        for (auto iter = fineLayout.begin(); iter.ok(); ++iter)
        {
            for (int dir = 0; dir < DIM; dir++)
            {
                BoxData<double> flux(iter.box());
                FR.incrementFine(flux, *iter, 1, dir);
            }
        }
        for (auto iter = crseLayout.begin(); iter.ok(); ++iter)
        {
            for (int dir = 0; dir < DIM; dir++)
            {
                BoxData<double> flux(iter.box());
                FR.incrementCoarse(flux, *iter, 1, dir);
            }
        }

        //pout() << "Error: " << err[nn] << std::endl;
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

