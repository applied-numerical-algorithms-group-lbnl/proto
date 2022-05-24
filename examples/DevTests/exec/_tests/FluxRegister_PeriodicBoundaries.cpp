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
    int testNum = 0;

    args.add("domainSize",      domainSize);
    args.add("physDomainSize",  physDomainSize);
    args.add("boxSize",         boxSize);
    args.add("refRatio",        refRatio);
    args.add("testNum",         testNum);
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

    if (procID() == 0) { std::cout << "dx = " << dx[0] << std::endl; }
    crseLayout.print();
    fineLayout.print();

    if (testNum == 0)
    {
        if (procID() == 0)
        {
            std::cout << "TEST 0: MEMTYPE = HOST" << std::endl;
        }
    
        LevelFluxRegister<double, 1, HOST> FR_CF(crseLayout, fineLayout, Point::Ones(refRatio), dx);
        LevelFluxRegister<double, 1, HOST> FR_C(crseLayout, fineLayout, Point::Ones(refRatio), dx);
        LevelFluxRegister<double, 1, HOST> FR_F(crseLayout, fineLayout, Point::Ones(refRatio), dx);
        for (auto iter = fineLayout.begin(); iter.ok(); ++iter)
        {
            for (int dir = 0; dir < DIM; dir++)
            {
                BoxData<double, 1, HOST> flux(iter.box().grow(dir, Side::Hi, 1));
                flux.setVal((*iter).global()+1);
                FR_F.incrementFine(flux, *iter, dx[dir], dir);
                FR_CF.incrementFine(flux, *iter, dx[dir], dir);
            }
        }
        for (auto iter = crseLayout.begin(); iter.ok(); ++iter)
        {
            for (int dir = 0; dir < DIM; dir++)
            {
                BoxData<double, 1, HOST> flux(iter.box().grow(dir, Side::Hi, 1));
                flux.setVal(((*iter).global()+1)*100);
                FR_C.incrementCoarse(flux, *iter, dx[dir], dir);
                FR_CF.incrementCoarse(flux, *iter, dx[dir], dir);
            }
        }

        h5.writeFluxRegister(dx, FR_CF, "FR_CF_0");
        h5.writeFluxRegister(dx, FR_C, "FR_C_0");
        h5.writeFluxRegister(dx, FR_F, "FR_F_0");
    } else if (testNum == 1)
    {
        if (procID() == 0)
        {
            std::cout << "TEST 1: MEMTYPE = DEVICE" << std::endl;
        }
    
        LevelFluxRegister<double, 1, DEVICE> FR_CF(crseLayout, fineLayout, Point::Ones(refRatio), dx);
        LevelFluxRegister<double, 1, DEVICE> FR_C(crseLayout, fineLayout, Point::Ones(refRatio), dx);
        LevelFluxRegister<double, 1, DEVICE> FR_F(crseLayout, fineLayout, Point::Ones(refRatio), dx);
        for (auto iter = fineLayout.begin(); iter.ok(); ++iter)
        {
            for (int dir = 0; dir < DIM; dir++)
            {
                BoxData<double, 1, DEVICE> flux(iter.box().grow(dir, Side::Hi, 1));
                flux.setVal((*iter).global()+1);
                FR_F.incrementFine(flux, *iter, dx[dir], dir);
                FR_CF.incrementFine(flux, *iter, dx[dir], dir);
            }
        }
        for (auto iter = crseLayout.begin(); iter.ok(); ++iter)
        {
            for (int dir = 0; dir < DIM; dir++)
            {
                BoxData<double, 1, DEVICE> flux(iter.box().grow(dir, Side::Hi, 1));
                flux.setVal(((*iter).global()+1)*100);
                FR_C.incrementCoarse(flux, *iter, dx[dir], dir);
                FR_CF.incrementCoarse(flux, *iter, dx[dir], dir);
            }
        }

        h5.writeFluxRegister(dx, FR_CF, "FR_CF_1");
        h5.writeFluxRegister(dx, FR_C, "FR_C_1");
        h5.writeFluxRegister(dx, FR_F, "FR_F_1");
    }

    #ifdef PR_MPI
    MPI_Finalize();
    #endif
    return 0;
}

