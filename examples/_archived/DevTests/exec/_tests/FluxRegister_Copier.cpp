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
    double physDomainSize = 1.0;
    int boxSize = 16;
    int refRatio = 2;
    std::array<bool, DIM> periodicity;
    periodicity.fill(true);

    InputArgs args;
    args.add("domainSize",      domainSize);
    args.add("physDomainSize",  physDomainSize);
    args.add("boxSize",         boxSize);
    args.add("refRatio",        refRatio);
    args.add("periodic_x",      periodicity[0]);
    args.add("periodic_y",      periodicity[1]);
#if DIM > 2
    args.add("periodic_z",      periodicity[2]);
#endif
    args.parse(argc, argv);
    args.print();

    // BUILD GRIDS
    std::array<double, DIM> dx;
    dx.fill(physDomainSize / domainSize);

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
    LevelFluxRegister<double> fluxRegister(crseLayout, fineLayout, Point::Ones(refRatio));
    h5.writeFluxRegister(dx, fluxRegister, "FLUX_0");
    for (auto fiter = fineLayout.begin(); fiter.ok(); ++fiter)
    {
        for (int dir = 0; dir < DIM; dir++)
        {
            BoxData<double> flux(fiter.box().grow(dir, Side::Hi, 1));
            flux.setVal(17);
            fluxRegister.incrementFine(flux, *fiter, 1.0, dir);
        }
    }
    h5.writeFluxRegister(dx, fluxRegister, "FLUX_1");
    for (auto citer = crseLayout.begin(); citer.ok(); ++citer)
    {
        for (int dir = 0; dir < DIM; dir++)
        {
            BoxData<double> flux(citer.box().grow(dir, Side::Hi, 1));
            flux.setVal(5);
            fluxRegister.incrementCoarse(flux, *citer, 1.0, dir);
        }
    }
    h5.writeFluxRegister(dx, fluxRegister, "FLUX_2");

    LevelBoxData<double> data(crseLayout, Point::Ones());
    data.setToZero();
    fluxRegister.reflux(data, 1.0);
    h5.writeLevel(dx, data, "DATA"); 

    #ifdef PR_MPI
    MPI_Finalize();
    #endif
    return 0;
}

