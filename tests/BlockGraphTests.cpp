#include <gtest/gtest.h>
#include "ProtoMMB.H"

using namespace Proto;

TEST(BlockGraph, BasicFunctions) {
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
    for (int ii = 0; ii < numBlocks; ii++)
    {
        graph.addBoundary(ii, (ii+1) % numBlocks,
            Point::Basis(0), Side::Hi,
            CoordPermutation::ccw());
    }
    Box K = Box::Kernel(3);
    for (int ii = 0; ii < numBlocks; ii++)
    {
        const Block& block = graph[ii];
        for (auto dir : K)
        {
            switch (dir)
            {
                case Point::Basis(0):
                case Point::Basis(0) + Point::Basis(1, -1):
                    EXPECT_EQ(block.numBoundaries(dir), 1);
                    EXPECT_EQ(block[dir][0]->id(), (ii+1)%numBlocks);
                    break;
                case Point::Basis(1):
                case Point::Basis(0,-1) + Point::Basis(1):
                    EXPECT_EQ(block.numBoundaries(dir), 1);
                    EXPECT_EQ(block[dir][0]->id(), (ii+numBlocks-1) % numBlocks);
                    break;
                case Point::Basis(0) + Point::Basis(1):
                    EXPECT_EQ(block.numBoundaries(dir), 2);
                    break;
                default:
                    EXPECT_EQ(block.numBoundaries(dir), 0);
                    break;
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
