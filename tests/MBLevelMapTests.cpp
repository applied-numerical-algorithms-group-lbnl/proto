#include <gtest/gtest.h>
#include "Proto.H"
#include "Lambdas.H"
#include "MBLevelMap_Shear.H"
#include "MBLevelMap_CubeSphereShell.H"

using namespace Proto;

TEST(MBMapTests, ShearMap) {
    int domainSize = 8;
    int boxSize = 8;
    HDF5Handler h5;

    auto domain = buildShear(domainSize);
    Point boxSizeVect = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);

    Array<Point, DIM+1> ghost;
    ghost.fill(Point::Ones(4));
    ghost[0] = Point::Ones(1);

    // initialize map
    MBLevelMap_Shear<HOST> map;
    map.define(layout, ghost);

    
#if PR_VERBOSE > 0
    h5.writeMBLevel({"X", "Y", "Z"}, map, map.map(), "MBMapTests_ShearMap_X");
#endif
}
#if DIM > 2
TEST(MBMapTests, CubeSphereShell) {
    int domainSize = 8;
    int boxSize = 8;
    int thickness = 8;
    int radialDir = 0;
    HDF5Handler h5;

    auto domain = buildCubeSphereShell(domainSize, thickness, radialDir);
    Point boxSizeVect = Point::Ones(boxSize);
    boxSizeVect[radialDir] = thickness;
    MBDisjointBoxLayout layout(domain, boxSizeVect);

    Array<Point, DIM+1> ghost;
    ghost.fill(Point::Ones(4));
    ghost[0] = Point::Ones(1);

    // initialize map
    MBLevelMap_CubeSphereShell<HOST> map;
    map.define(layout, ghost);
    
#if PR_VERBOSE > 0
    h5.writeMBLevel({"X", "Y", "Z"}, map, map.map(), "MBMapTests_CubeSphereMap_X");
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
