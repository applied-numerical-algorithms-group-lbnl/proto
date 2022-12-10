#include <gtest/gtest.h>
#include "Proto.H"

#define XPOINT_SIZE 5


using namespace Proto;
MBGraph buildXPoint()
{
    MBGraph graph(XPOINT_SIZE);
    auto CCW = CoordPermutation::ccw();
    for (int ii = 0; ii < XPOINT_SIZE; ii++)
    {
        graph.defineBoundary(ii, (ii+1) % XPOINT_SIZE, 0, Side::Hi, CCW);
    }
    graph.fixRotations();
    return graph;
}

TEST(MBGraph, XPointConnectivity) {
    int numBlocks = XPOINT_SIZE;
    MBGraph graph = buildXPoint();
    for (int bi = 0; bi < numBlocks; bi++)
    {
        for (int bj = 0; bj < numBlocks; bj++)
        {
            Point dir_ij = graph.connectivity(bi, bj);
            Point dir_ji = graph.connectivity(bj, bi);
            int diff = bj - bi;
            if (diff == 0)
            {
                EXPECT_EQ(dir_ij, Point::Zeros());
                EXPECT_EQ(dir_ji, Point::Zeros());
            }
            else if (diff == 1 || diff == -(XPOINT_SIZE-1))
            {
                EXPECT_EQ(dir_ij, Point::Basis(0));
                EXPECT_EQ(dir_ji, Point::Basis(1));
            } else if (diff == -1 || diff == XPOINT_SIZE-1) {
                EXPECT_EQ(dir_ij, Point::Basis(1));
                EXPECT_EQ(dir_ji, Point::Basis(0));
            } else {
                EXPECT_EQ(dir_ij, Point::Basis(0) + Point::Basis(1));
                EXPECT_EQ(dir_ji, Point::Basis(0) + Point::Basis(1));
            }
        }
    }
}

TEST(MBGraph, XPointReverseArc) {
    int numBlocks = XPOINT_SIZE;
    MBGraph graph = buildXPoint();
    Point x = Point::Basis(0);
    Point y = Point::Basis(1);
    Box K = Box::Kernel(1);
    for (int bi = 0; bi < numBlocks; bi++)
    {
        for (int bj = 0; bj < numBlocks; bj++)
        {
            int diff = bj - bi;
            if (diff == 0) { continue; }
            else if (diff == 1 || diff == -(XPOINT_SIZE-1))
            {
                Box orth = K.grow(-x);
                auto R = graph.rotation(bi, x);
                for (auto dir : orth)
                {
                    Point arc = dir + x;
                    Point revArc = R(dir - x);
                    EXPECT_EQ(graph.reverseArc(bi,  bj, arc), revArc);
                }
            } else if (diff == -1 || diff == XPOINT_SIZE-1) {
                Box orth = K.grow(-y);
                auto R = graph.rotation(bi, y);
                for (auto dir : orth)
                {
                    Point arc = dir + y;
                    Point revArc = R(dir - y);
                    EXPECT_EQ(graph.reverseArc(bi,  bj, arc), revArc);
                }
            } else {
                
            }
        }
    }
}

TEST(MBGraph, XPointAdjacent) {
    int numBlocks = XPOINT_SIZE;
    MBGraph graph = buildXPoint();
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

TEST(MBGraph, XPointFullConnectivity) {
    int numBlocks = XPOINT_SIZE;
    MBGraph graph = buildXPoint();
    for (int src = 0; src < numBlocks; src++)
    {
        for (int dst = 0; dst < numBlocks; dst++)
        {
            auto connectivity = graph.fullConnectivity(src, dst);
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

TEST(MBGraph, XPointBoundaries) {
    int numBlocks = XPOINT_SIZE;
    MBGraph graph = buildXPoint();
    for (int bi = 0; bi < numBlocks; bi++)
    {
        Box K = Box::Kernel(1);
        for (auto dir : K)
        {
            if (dir == Point::Zeros()) { continue; }
            auto srcBounds = graph.boundaries(bi, dir);
            if (dir == Point::Basis(0) || dir == Point::Basis(1))
            {
                EXPECT_EQ(srcBounds.size(), 1);
            } else if (dir == (Point::Basis(0) + Point::Basis(1)))
            {
                EXPECT_EQ(srcBounds.size(), numBlocks-3);
            } else {
                EXPECT_EQ(srcBounds.size(), 0);
            }
            for (auto srcArc : srcBounds)
            {
                EXPECT_EQ(bi, srcArc.srcBlock);
                EXPECT_EQ(dir, srcArc.srcToDst);
                unsigned int dstBlock = srcArc.dstBlock;
                auto dstBounds = graph.boundaries(dstBlock, srcArc.dstToSrc);
                bool foundReverseArc = false;
                for (auto dstArc : dstBounds)
                {
                    if (dstArc.dstBlock == bi)
                    {
                        EXPECT_FALSE(foundReverseArc);
                        EXPECT_EQ(srcArc.dstToSrc, dstArc.srcToDst);
                        EXPECT_EQ(srcArc.srcToDst, dstArc.dstToSrc);
                        foundReverseArc = true;
                        break;
                    }
                }
            }
        }
    }
}

TEST(MBGraph, CubedSphere) {
    int numBlocks = 2*DIM+1;
    MBGraph graph(numBlocks);
    CoordPermutation I;
    int id = 1;
    for (int di = 0; di < DIM; di++)
    {
        SideIterator si;
        for (si.begin(); si.ok(); ++si)
        {
            graph.defineBoundary(0, id, di, *si, I);
            id++;
        }
    }
    for (int di = 0; di < DIM; di++)
    {
        for (int dj = di+1; dj < DIM; dj++)
        {
            CoordPermutation R = CoordPermutation::rotatePlane(di, dj);
            Point dir = Point::Basis(di);
            for (int ii = 0; ii < 4; ii++)
            {
                Point adjDir = R(dir);
                auto srcBlock = graph.adjacent(0, dir);
                auto dstBlock = graph.adjacent(0, adjDir);
                if (graph.adjacent(srcBlock, adjDir) == graph.size())
                {
                    graph.defineBoundary(srcBlock, dstBlock, adjDir, R);
                }
                dir = adjDir;
            }
        }
    }
    graph.fixRotations();
    auto dirs = codimDirs(1);
    for (auto di : dirs)
    {
        auto bi = graph.adjacent(0, di);
        auto conn_out = graph.fullConnectivity(0, bi);
        auto conn_in  = graph.fullConnectivity(bi, 0);
        EXPECT_EQ(conn_out.size(), pow(3,DIM-1));
        EXPECT_EQ(conn_in.size(), pow(3,DIM-1));
        for (auto dj : dirs)
        {
            if (di.dot(dj) == 0)
            {
                auto bj = graph.adjacent(0, dj);
                auto conn_ij = graph.fullConnectivity(bi, bj);
                EXPECT_EQ(conn_ij.size(), pow(3,DIM-1));
            }
        }
    }
}
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
