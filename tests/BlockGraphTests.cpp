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

TEST(BlockGraph, PlanarXPoint) {
    int numBlocks = XPOINT_SIZE;
    BlockGraph graph = buildXPoint();
    Box K = Box::Kernel(1);
    for (int ii = 0; ii < numBlocks; ii++)
    {
        const BlockGraphNode& block = graph[ii];
        for (auto dir : K)
        {
#if DIM == 2
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
#elif DIM == 3
            if (
                    dir == Point::Basis(0)
                 || dir == Point::Basis(0) + Point::Basis(2,+1)
                 || dir == Point::Basis(0) + Point::Basis(2,-1)
                 || dir == Point::Basis(0) + Point::Basis(1,-1)
                 || dir == Point::Basis(0) + Point::Basis(1,-1) + Point::Basis(2,+1)
                 || dir == Point::Basis(0) + Point::Basis(1,-1) + Point::Basis(2,-1))
            {
                EXPECT_EQ(block.numBoundaries(dir), 1);
                EXPECT_EQ((*block.boundaries(dir).begin())->id(), (ii+1)%numBlocks);
            } else if(
                    dir == Point::Basis(1)
                 || dir == Point::Basis(1) + Point::Basis(2,+1)
                 || dir == Point::Basis(1) + Point::Basis(2,-1)
                 || dir == Point::Basis(1) + Point::Basis(0,-1)
                 || dir == Point::Basis(1) + Point::Basis(0,-1) + Point::Basis(2,+1)
                 || dir == Point::Basis(1) + Point::Basis(0,-1) + Point::Basis(2,-1))
            {
                EXPECT_EQ(block.numBoundaries(dir), 1);
                EXPECT_EQ((*block.boundaries(dir).begin())->id(), (ii+numBlocks-1) % numBlocks);
            } else if (
                    dir == Point::Basis(0) + Point::Basis(1)
                 || dir == Point::Basis(0) + Point::Basis(1) + Point::Basis(2,1)
                 || dir == Point::Basis(0) + Point::Basis(1) + Point::Basis(2,-1))
            {
                EXPECT_EQ(block.numBoundaries(dir), 2);
            } else {
                EXPECT_EQ(block.numBoundaries(dir), 0);
            }
#else
            EXPECT_TRUE(false)
#endif
        }
    }
}

TEST(BlockGraph, Adjacent) {
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

TEST(BlockGraph, Connectivity) {
    int numBlocks = XPOINT_SIZE;
    BlockGraph graph = buildXPoint();
    for (int src = 0; src < numBlocks; src++)
    {
        for (int dst = 0; dst < numBlocks; dst++)
        {
            auto connectivity = graph.connectivity(src, dst);
            std::cout << "Checking connectivity from " << src << " to " << dst << std::endl;
            for (auto item : connectivity)
            {
                std::cout << item << ", ";
            }
            std::cout << std::endl;
            Box dirs;
            if (src == dst)
            {
                // do nothing
            } else if (dst - src == 1) {
                dirs = Box::Kernel(1).grow(0, -1).shift(Point::Basis(0));
            } else if (src - dst == 1) {
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
