#include "ProtoAMR.H"
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

    int domainSize = 64;
    double physDomainSize = 1.0;
    int boxSize = 16;
    int refRatio = 2;
    std::array<bool, DIM> periodicity;
    periodicity.fill(true);

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

    std::array<double, DIM> dx;
    dx.fill(physDomainSize / domainSize);

    Point boxSizeV = Point::Ones(boxSize);
    Box domainBox = Box::Cube(domainSize);
    ProblemDomain crseDomain(domainBox, periodicity);
    DisjointBoxLayout crseLayout(crseDomain, boxSizeV);

    Box fineDomainBox = Box::Cube(domainSize/2).refine(refRatio);
    ProblemDomain fineDomain = crseDomain.refine(refRatio);
    DisjointBoxLayout fineLayout(fineDomain, fineDomainBox, boxSizeV);

    crseLayout.print();
    fineLayout.print();

    LevelFluxRegister<double> FR(crseLayout, fineLayout, Point::Ones(refRatio), dx);
    FR.print();
    h5.writeFluxRegister(dx, FR, "FR_0");
    for (auto iter = fineLayout.begin(); iter.ok(); ++iter)
    {
        for (int dir = 0; dir < DIM; dir++)
        {
            BoxData<double> flux(iter.box());
            flux.setVal(1);
            FR.incrementFine(flux, *iter, 1, dir);
        }
    }
    h5.writeFluxRegister(dx, FR, "FR_1");
    for (auto iter = crseLayout.begin(); iter.ok(); ++iter)
    {
        for (int dir = 0; dir < DIM; dir++)
        {
            BoxData<double> flux(iter.box());
            flux.setVal(10);
            FR.incrementCoarse(flux, *iter, 1, dir);
        }
    }
       
    h5.writeFluxRegister(dx, FR, "FR_2");
        
    #ifdef PR_MPI
    MPI_Finalize();
    #endif
    return 0;
}

