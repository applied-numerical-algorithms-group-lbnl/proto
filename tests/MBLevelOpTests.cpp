#include <gtest/gtest.h>
#include "Proto.H"
#include "Lambdas.H"
#include "BoxOp_TestMBLaplace.H"

using namespace Proto;

TEST(MBLevelOp, Laplace) {
    
    HDF5Handler h5;
    int domainSize = 64;
    int boxSize = 32;
    auto domain = buildXPoint(domainSize);
    Point boxSizeVect = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);
    Array<Point, DIM+1> ghost;
    ghost.fill(Point::Zeros());
    ghost[0] = Point::Ones(3);
    Point boundGhost = Point::Ones();
   
    // Map
    MBMap<XPointMapRigid_t> map(XPointMapRigid, layout, ghost, boundGhost);
    
    // Laplace Operator
    MBLevelOp<BoxOp_TestMBLaplace, double, XPointMapRigid_t> op(map);
   
#if 0
    // Interpolation Operator
    std::vector<Point> footprint;
    for (auto pi : Box::Kernel(1))
    {
        footprint.push_back(pi);
    }
    for (int dir = 0; dir < DIM; dir++)
    {
        footprint.push_back(Point::Basis(dir,2));
        footprint.push_back(Point::Basis(dir,-2));
    }
    MBInterpOp interp(map, footprint, ghost[0], 4);
#endif
    // Initialize Data
    MBLevelBoxData<double, 1, HOST> hostSrc(layout, ghost);
    MBLevelBoxData<double, 1, HOST> hostDst(layout, ghost);

    Point k{1,2,3,4,5,6};
    Array<double, DIM> offset{1,1,1,1,1,1};
    hostSrc.initialize(f_phiM, map, k, offset);
#if 0
    hostSrc.fillBoundaries();
    interp.apply(hostSrc, hostSrc);
#endif

    // Apply Operator
    op(hostDst, hostSrc);
    h5.writeMBLevel({"phi"}, map, hostSrc, "PHI");
    h5.writeMBLevel({"Lphi"}, map, hostDst, "LPHI");
}

int main(int argc, char *argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
#ifdef PR_MPI
    MPI_Init(&argc, &argv);
#endif
    int result = RUN_ALL_TESTS();
#ifdef PR_MPI
    MPI_Finalize();
#endif
    return result;
}
