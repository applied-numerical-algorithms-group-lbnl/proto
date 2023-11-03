#include <gtest/gtest.h>
#include "Proto.H"
#include "Lambdas.H"
#include "BoxOp_MBLaplace.H"
#include "MBMap_XPointRigid.H"

using namespace Proto;

TEST(MBMultigridTests, Construction) {
    int domainSize = 32;
    int boxSize = 16;
    int numBlocks = 5;
    int numLevels = 3;

    auto domain = buildXPoint(domainSize, numBlocks);
    MBDisjointBoxLayout layout(domain, Point::Ones(boxSize));
    std::vector<MBLevelBoxData<double, 1, HOST>> v;

    Array<Point,DIM+1> ghost;
    ghost.fill(Point::Zeros());
    MBMultigrid<BoxOp_MBLaplace, MBMap_XPointRigid, double> mg(layout, Point::Ones(2), numLevels); 
}

int main(int argc, char *argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
#ifdef PR_MPI
    MPI_Init(&argc, &argv);
#endif
    PR_TIMER_SETFILE("MBMultigridTests_DIM" + to_string(DIM) + "_NProc" + to_string(numProc())
            + ".time.table");
    int result = RUN_ALL_TESTS();
    PR_TIMER_REPORT();
#ifdef PR_MPI
    MPI_Finalize();
#endif
    return result;
}
