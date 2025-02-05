#include <gtest/gtest.h>
#include "Proto.H"
#include "TestFunctions.H"
using namespace Proto;

TEST(Suite, Test) {
    int domainSize = 16;
    int boxSize = 16;
    int numBlocks = 3;
    int ghostSize = 2;
    auto domain = buildXPoint(domainSize, numBlocks);
    Point boxSizeVect = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, Point::Ones(boxSize));

    MBLevelMap<MBMap_XPointRigid, HOST> map;
    map.define(layout, Point::Ones(ghostSize));
    for (auto bi = 0; bi < numBlocks; bi++)
    {
        map[bi].setNumBlocks(numBlocks);
    }
    map.initialize();

    MBInterpOp interp(map);
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
