#include <gtest/gtest.h>
#include "ProtoMMB.H"

#define XPOINT_SIZE 5

using namespace Proto;

BlockGraph buildXPoint()
{
    BlockGraph graph(XPOINT_SIZE);
    auto CCW = CoordPermutation::ccw();
    for (int ii = 0; ii < XPOINT_SIZE; ii++)
    {
        graph.addBoundary(ii, (ii+1) % XPOINT_SIZE, 0, Side::Hi, CCW);
    }
    graph.print();
    return graph;
}

TEST(BlockGraph, XPointAdjacent) {
    int numBlocks = XPOINT_SIZE;
    BlockGraph graph = buildXPoint();
    for (int ii = 0; ii < numBlocks; ii++)
    {
        for (int dir = 0; dir < DIM; dir++)
        {
            SideIterator siter;
            for (siter.begin(); siter.ok(); ++siter)
            {
                if (dir == 0 && *siter == Side::Hi)
                {
                    EXPECT_EQ(graph.adjacent(ii, dir, *siter), (ii+1) % numBlocks);
                } else if (dir == 1 && *siter == Side::Hi)
                {
                    EXPECT_EQ(graph.adjacent(ii, dir, *siter), (ii + numBlocks - 1) % numBlocks);
                } else {
                    EXPECT_EQ(graph.adjacent(ii, dir, *siter), numBlocks);
                }
            }
        }
    }
}

TEST(BlockGraph, XPointConnectivity) {
    int numBlocks = XPOINT_SIZE;
    BlockGraph graph = buildXPoint();
    for (int src = 0; src < numBlocks; src++)
    {
        for (int dst = 0; dst < numBlocks; dst++)
        {
            auto connectivity = graph.connectivity(src, dst);
            Box dirs;
            if (dst == src)
            {
                // do nothing
            } else if (dst == (src+1)%numBlocks) {
                dirs = Box::Kernel(1).grow(0, -1).shift(Point::Basis(0));
            } else if (dst == (src+numBlocks-1)%numBlocks) {
                dirs = Box::Kernel(1).grow(1, -1).shift(Point::Basis(1));
            } else {
                dirs = Box::Kernel(1).edge(Point::Basis(0) + Point::Basis(1));
            }
            EXPECT_EQ(connectivity.size(), dirs.size());
            for (auto dir : dirs)
            {
                EXPECT_NE(connectivity.find(dir), connectivity.end());
            }
        }
    }
}
/*
TEST(BlockGraph, CubedSphere) {
    int numBlocks = 2*DIM+1;
    BlockGraph graph(numBlocks);
    CoordPermutation I;
    int id = 1;
    for (int dir = 0; dir < DIM; dir++)
    {
        SideIterator siter;
        for (siter.begin(); siter.ok(); ++siter)
        {
            graph.addBoundary(0, id, dir, *siter, I);
        }
    }
    for (auto dir : codimDirs(2))
    {
        graph.closeCircuit(0,dir);
    }
}
*/
int main(int argc, char *argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
#ifdef PR_MPI
    MPI_Init(&argc, &argv);
#endif
    std::cout << "DIM = " << DIM << std::endl;
    int result = RUN_ALL_TESTS();
#ifdef PR_MPI
    MPI_Finalize();
#endif
    return result;
}
