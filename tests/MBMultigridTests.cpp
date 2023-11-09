#include <gtest/gtest.h>
#include "Proto.H"
#include "Lambdas.H"
#include "BoxOp_MBLaplace.H"
#include "MBMap_XPointRigid.H"

using namespace Proto;

#if 1
TEST(MBMultigridTests, Construction) {
    typedef BoxOp_MBLaplace<double> OP;
    int domainSize = 32;
    int boxSize = 16;
    int numBlocks = 5;
    int numLevels = 3;
    Point refRatio = Point::Ones(2);
    std::vector<Point> refRatios;
    for (int bi = 0; bi < numBlocks; bi++)
    {
        refRatios.push_back(refRatio);
    }
    auto domain = buildXPoint(domainSize, numBlocks);
    MBDisjointBoxLayout layout(domain, Point::Ones(boxSize));

    MBMultigrid<BoxOp_MBLaplace, MBMap_XPointRigid, double> mg(layout, refRatio, numLevels); 
    EXPECT_TRUE(mg.validate(layout, refRatios, numLevels));


    Array<Point,DIM+1> ghost;
    Array<Point,DIM+1> noGhost;
    ghost.fill(OP::ghost());
    noGhost.fill(Point::Zeros());
    MBLevelBoxData<double, 1, HOST> phi(layout, ghost);
    MBLevelBoxData<double, 1, HOST> rhs(layout, noGhost);

    mg.solve(phi, rhs, 10, 1e-10);
}
#endif

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
