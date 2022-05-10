#include <gtest/gtest.h>
#include "ProtoMMB.H"

using namespace Proto;

TEST(BlockGraph, BasicFunctions) {
    CoordPermutation R;
    BlockGraphNode block(0);
    int numBlocks = 3;
    BlockGraph graph(numBlocks);
    EXPECT_EQ(graph.size(), numBlocks);
    for (int ii = 0; ii < numBlocks; ii++)
    {
        EXPECT_EQ(graph[ii].id(), ii);
    }
}
TEST(BlockGraph, SquareCircle) {
    int numBlocks = 3;
    BlockGraph graph(numBlocks);
    auto CCW = CoordPermutation::ccw();
    for (int ii = 0; ii < numBlocks; ii++)
    {
        graph.addBoundary(ii, (ii+1) % numBlocks, 0, Side::Hi, CCW);
    }
    Box K = Box::Kernel(3);
    for (int ii = 0; ii < numBlocks; ii++)
    {
        const BlockGraphNode& block = graph[ii];
        for (auto dir : K)
        {
            if (dir == Point::Basis(0) || dir == Point::Basis(0) + Point::Basis(1,-1))
            {
                EXPECT_EQ(block.numBoundaries(dir), 1);
                EXPECT_EQ(block.boundaries(dir)[0]->id(), (ii+1)%numBlocks);
            } else if(dir == Point::Basis(1) || dir == Point::Basis(0,-1) + Point::Basis(1))
            {
                EXPECT_EQ(block.numBoundaries(dir), 1);
                EXPECT_EQ(block.boundaries(dir)[0]->id(), (ii+numBlocks-1) % numBlocks);
            } else if (dir == Point::Basis(0) + Point::Basis(1))
            {
                EXPECT_EQ(block.numBoundaries(dir), 2);
            } else {
                EXPECT_EQ(block.numBoundaries(dir), 0);
            }
        }
    }
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
