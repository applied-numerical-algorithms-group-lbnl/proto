#include <gtest/gtest.h>
#include "Proto.H"
#include "TestFunctions.H"
#include "MBMap_Shear.H"
#include "BoxOp_MBLaplace.H"

using namespace Proto;

TEST(MBLevelRK4Tests, Construction) {
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
    MBLevelMap<MBMap_Shear, HOST> map;
    map.define(layout, ghost);

    MBInterpOp interp(map, 4);
    MBLevelRK4<BoxOp_MBLaplace, MBMap_Shear, double> rk4(map, interp);
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
