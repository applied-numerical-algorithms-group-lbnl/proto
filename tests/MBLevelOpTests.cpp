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
    int nghost = 4;
    ghost.fill(Point::Ones(nghost+2));
    ghost[0] = Point::Ones(nghost);
    Point boundGhost = Point::Ones();
   
    // Map
    MBMap<XPointMapRigid_t> map(XPointMapRigid, layout, ghost, boundGhost);
    
    // Laplace Operator
    MBLevelOp<BoxOp_TestMBLaplace, double, XPointMapRigid_t> op(map);
   
#if 1
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
    MBLevelBoxData<double, 1, HOST> hostDstNorm(layout, ghost);

    Point k{1,1,1,1,1,1};
    Array<double, DIM> offset{1,1,1,1,1,1};
    hostSrc.initialize(f_phiM, map, k, offset);
#if 1
    hostSrc.fillBoundaries();
    interp.apply(hostSrc, hostSrc);
#endif

    // Apply Operator
    op(hostDst, hostSrc);
#if PR_VERBOSE > 0
    h5.writeMBLevel({"phi"}, map, hostSrc, "MBLevelOpTests_Phi_0");
    h5.writeMBLevel({"Lphi"}, map, hostDst, "MBLevelOpTests_LPhi_0");
#endif
    hostDst.fillBoundaries();
    interp.apply(hostDst, hostDst);
#if PR_VERBOSE > 0
    h5.writeMBLevel({"phi"}, map, hostSrc, "MBLevelOpTests_Phi_1");
    h5.writeMBLevel({"Lphi"}, map, hostDst, "MBLevelOpTests_LPhi_1");
#endif
    for (auto iter : layout)
    {
        auto& hostDstNorm_i = hostDstNorm[iter];
        auto& hostDst_i = hostDst[iter];
        const auto& J_i = map.jacobian()[iter];
        Operator::cellQuotient(hostDstNorm_i, hostDst_i, J_i);
    }
#if PR_VERBOSE > 0
    h5.writeMBLevel({"LphiNorm"}, map, hostDstNorm, "MBLevelOpTests_LPhiNorm");
#endif
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
