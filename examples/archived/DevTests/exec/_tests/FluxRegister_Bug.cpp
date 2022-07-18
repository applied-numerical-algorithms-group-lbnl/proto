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

    int domainSize = 32;
    double physDomainSize = domainSize;
    int boxSize = 32;
    int refRatio = 4;
    std::array<bool, DIM> periodicity;
    for (int dir = 0; dir < DIM; dir++) { periodicity[dir] = true; }

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
    double cdx = physDomainSize / domainSize;
    double fdx = cdx / refRatio;
    
    std::array<double, DIM> cdxVect;
    cdxVect.fill(cdx);

    // coarse layout
    Point boxSizeV = Point::Ones(boxSize);
    Box domainBox = Box::Cube(domainSize);
    ProblemDomain crseDomain(domainBox, periodicity);
    DisjointBoxLayout crseLayout(crseDomain, boxSizeV);

    // fine layout
    std::vector<Point> finePatches;
    Box b = Box::Cube(3).shift(Point::Ones());
    for (auto biter = b.begin(); biter != b.end(); ++biter)
    {
        if (*biter == Point::Ones(3)){continue;}
        finePatches.push_back(*biter);
    }

    ProblemDomain fineDomain = crseDomain.refine(Point::Ones(refRatio));
    DisjointBoxLayout fineLayout(fineDomain, finePatches, boxSizeV);
    
    LevelBoxData<double> fineData(fineLayout, Point::Zeros());
    LevelBoxData<double> crseData(crseLayout, Point::Ones());

    crseLayout.print();
    fineLayout.print();

    h5.writeLevel(fdx, fineData, "FINE_DATA");
    h5.writeLevel(cdx, crseData, "CRSE_DATA");
    
    LevelFluxRegister<double> fineFlux(crseLayout, fineLayout, Point::Ones(refRatio));

    for (auto iter = fineLayout.begin(); iter.ok(); ++iter)
    {
        int index = (*iter).global();
        for (int dir = 0; dir < DIM; dir++)
        {
            Box fluxBox = iter.box().grow(dir, Side::Hi, 1);
            BoxData<double> flux(fluxBox);
            flux.setVal(1000 + 100*index + 10*dir + procID());
            fineFlux.incrementFine(flux, *iter, 1.0, dir);
        }
    }
    h5.writeFluxRegister(cdxVect, fineFlux, "FINE_FLUX");

    #ifdef PR_MPI
    MPI_Finalize();
    #endif
    return 0;
}

