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

    int domainSize = 64;
    int boxSize = 16;
    std::array<bool, DIM> periodicity;
    for (int dir = 0; dir < DIM; dir++) { periodicity[dir] = true; }

    args.set("periodic_x", &periodicity[0]);
    args.set("periodic_y", &periodicity[1]);

    Point boxSizeV = Point::Ones(boxSize);
    Box domainBox = Box::Cube(domainSize);
    ProblemDomain domain(domainBox, periodicity);
    DisjointBoxLayout layout(domain, boxSizeV);
    
    DisjointBoxLayout cfLayout = layout.coarsen(Point::Ones(2));

    Box crseDomainBox = Box::Cube(domainSize / 2);
    ProblemDomain crseDomain(crseDomainBox, periodicity);
    DisjointBoxLayout crseLayout(crseDomain, boxSizeV);
    
    DisjointBoxLayout crseLayoutAlt(crseDomain, boxSizeV / 2);

    std::vector<Point> patches_0;
    patches_0.push_back(Point::Zeros());
    std::vector<Point> patches_1;
    patches_1.push_back(Point::Ones());
    DisjointBoxLayout singlePatchLayout_0(domain, patches_0, boxSizeV);
    DisjointBoxLayout singlePatchLayout_1(domain, patches_1, boxSizeV);

    bool allTestsPass = true;
    bool compatible;

    pout() << "Checking Compatibility: " << std::endl;
    pout() << "===========================================================" << std::endl << std::endl;
    layout.print();
    cfLayout.print();
    compatible = layout.compatible(cfLayout); 
    pout() << "Are layouts compatible? " << compatible << std::endl;
    allTestsPass &= compatible; //should be true
    pout() << "===========================================================" << std::endl << std::endl;
    layout.print();
    crseLayout.print();
    compatible = layout.compatible(crseLayout);
    pout() << "Are layouts compatible? " << compatible << std::endl;
    allTestsPass &= (!compatible); //should be false
    pout() << "===========================================================" << std::endl << std::endl;
    layout.print();
    crseLayoutAlt.print();
    compatible = layout.compatible(crseLayoutAlt);
    pout() << "Are layouts compatible? " << compatible << std::endl;
    allTestsPass &= (compatible); //should be true
    pout() << "===========================================================" << std::endl << std::endl;
    singlePatchLayout_0.print();
    singlePatchLayout_1.print();
    compatible = singlePatchLayout_0.compatible(singlePatchLayout_1);
    pout() << "Are layouts compatible? " << compatible << std::endl;
    allTestsPass &= (!compatible); //should be false
    pout() << "===========================================================" << std::endl << std::endl;
    if (allTestsPass)
    {
        pout() << "ALL TESTS PASSED" << std::endl;
    } else {
        pout() << "SOME TESTS FAILED" << std::endl;
    }
     
    #ifdef PR_MPI
    MPI_Finalize();
    #endif
    return 0;
}

