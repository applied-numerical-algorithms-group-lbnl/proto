#include <gtest/gtest.h>
#include "Proto.H"
#include "Lambdas.H"
#include "BoxOp_TestMBLaplace.H"

using namespace Proto;

TEST(MBLevelOp, Laplace) {

    int domainSize = 8;
    int boxSize = 4;
    auto domain = buildXPoint(domainSize);
    Point boxSizeVect = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);
    Array<Point, DIM+1> ghost;
    ghost.fill(Point::Zeros());
    Point boundGhost = Point::Ones();
   
    MBMap<XPointMapRigid_t> map(XPointMapRigid, layout, ghost, boundGhost);

    MBLevelOp<BoxOp_TestMBLaplace, double, XPointMapRigid_t> op(map);
}

// TEST...
// TEST...
// TEST...
// TEST...

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
