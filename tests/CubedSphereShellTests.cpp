#include <gtest/gtest.h>
#include "Proto.H"
#include "Lambdas.H"

using namespace Proto;
#if DIM == 3
TEST(CubedSphereShell, BCInterpOp) {
#if PR_VERBOSE > 0
    HDF5Handler h5;
#endif
    int nSrcPhi = 360;
    int nSrcTheta = 180;

    int domainSize = 64;
    int boxSize = 32;
    int thickness = 1;
    int rCoord = CUBED_SPHERE_SHELL_RADIAL_COORD;
    int thetaCoord = (rCoord + 1) % 3;
    int phiCoord = (rCoord + 2) % 3;

    auto domain = CubedSphereShell::Domain(domainSize, thickness, rCoord);
    Point boxSizeVect = Point::Ones(boxSize);
    boxSizeVect[rCoord] = thickness;
    MBDisjointBoxLayout layout(domain, boxSizeVect);
    MBLevelBoxData<double, 1, HOST> dstData(layout, Point::Ones());
    dstData.setVal(0);

    // src data is a single block, single box MBLevelBoxData
    Point srcDomainSize;
    srcDomainSize[rCoord] = 1;
    srcDomainSize[phiCoord] = nSrcPhi;
    srcDomainSize[thetaCoord] = nSrcTheta;
    MBProblemDomain srcDomain(1);
    srcDomain.defineDomain(0,srcDomainSize);
    MBDisjointBoxLayout srcLayout(srcDomain, srcDomainSize);
    MBLevelBoxData<double, 1, HOST> srcData(srcLayout, Point::Zeros());
    srcData.setVal(1);

    auto map = CubedSphereShell::Map(layout, Point::Ones());

    std::vector<double> dtheta(nSrcTheta, M_PI/nSrcTheta);
    auto op = CubedSphereShell::BCInterpOp(map, srcLayout, dtheta, Side::Hi);
#if PR_VERBOSE > 0
    h5.writeMBLevel({"phi"}, srcData, "TestBCInterpOp_SrcData_0");
    h5.writeMBLevel({"phi"}, dstData, "TestBCInterpOp_DstData_0");
#endif

    op.apply(dstData, srcData);

#if PR_VERBOSE > 0
    h5.writeMBLevel({"phi"}, srcData, "TestBCInterpOp_SrcData_1");
    h5.writeMBLevel({"phi"}, dstData, "TestBCInterpOp_DstData_1");
#endif


}
#endif

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
