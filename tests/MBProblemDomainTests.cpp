#include <gtest/gtest.h>
#include "Proto.H"

#define XPOINT_SIZE 5

using namespace Proto;

namespace
{
    MBProblemDomain buildXPoint(int a_domainSize)
    {
        MBProblemDomain domain(XPOINT_SIZE);
        auto CCW = CoordPermutation::ccw();
        for (int ii = 0; ii < XPOINT_SIZE; ii++)
        {
            domain.defineBoundary(ii, (ii + 1) % XPOINT_SIZE, 0, Side::Hi, CCW);
        }
        for (int bi = 0; bi < XPOINT_SIZE; bi++)
        {
            domain.defineDomain(bi, Point::Ones(a_domainSize));
        }
        domain.close();
        return domain;
    }
}

TEST(MBProblemDomain, ConvertSimple)
{
    int boxSize = 4;
    for (auto dir : Box::Kernel(1))
    {
        if (dir.codim() != 1) { continue; }
        for (int axis = 0; axis < DIM; axis++)
        {
            if (dir[axis] != 0) { continue; }
            CoordPermutation R;
            for (int rot = 0; rot < 4; rot++)
            {
                MBProblemDomain domain(2);
                domain.defineBoundary(0,1,dir,R);
                domain.defineDomain(0,Point::Ones(boxSize));
                domain.defineDomain(1,Point::Ones(boxSize));
                const auto& graph = domain.graph();
                Point revDir = graph.reverseDir(0,1,dir);
                Point revArc = graph.reverseArc(0,1,dir);
                EXPECT_EQ(R(dir), -revDir);

                Box B0 = Box::Cube(boxSize).adjacent(dir);
                Box B1 = Box::Cube(boxSize).edge(revDir);
                Box B10 = domain.convertBox(B0,0,1);
                EXPECT_EQ(B10, B1);
                
                R = R*CoordPermutation::ccw(axis);
            }
        }
    }
}

TEST(MBProblemDomain, OnBlockBoundary)
{
    int domainSize = 4;
    auto domain = buildXPoint(domainSize);
    Box domainBox = Box::Cube(domainSize);
    std::vector<Box> blockBoundaries;
    std::vector<Box> domainBoundaries;

    for (auto dir : Box::Kernel(1))
    {
        if (dir == Point::Zeros())
        {
            continue; 
        }
        else if (dir == Point::X() || dir == Point::Y() || dir == Point::X() + Point::Y())
        {
            blockBoundaries.push_back(domainBox.edge(dir));
        } else {
            domainBoundaries.push_back(domainBox.edge(dir));
        }
    }

    BoxData<double, 1> blockBoundaryMask(domainBox);
    blockBoundaryMask.setVal(0);
    for (auto box : blockBoundaries)
    {
        blockBoundaryMask.setVal(1, box);
    }
    BoxData<double, 1> domainBoundaryMask(domainBox);
    domainBoundaryMask.setVal(0);
    for (auto box : domainBoundaries)
    {
        domainBoundaryMask.setVal(1, box);
    }

    for (BlockIndex block = 0; block < domain.numBlocks(); block++)
    {
        for (Point point : domainBox)
        {
            if (blockBoundaryMask(point) == 1)
            {
                EXPECT_TRUE(domain.isPointOnBlockBoundary(point, block));
            } else {
                EXPECT_FALSE(domain.isPointOnBlockBoundary(point, block));
            }
            if (domainBoundaryMask(point) == 1)
            {
                EXPECT_TRUE(domain.isPointOnDomainBoundary(point, block));
            } else {
                EXPECT_FALSE(domain.isPointOnDomainBoundary(point, block));
            }
        }
    }

}

TEST(MBProblemDomain, Convert) {
    int domainSize = 64;
    auto domain = buildXPoint(domainSize);
    Point x = Point::X();
    Point y = Point::Y();
    Point origin = Point::Zeros();
    Box domainBox = Box::Cube(domainSize);
    Box xAdj = domainBox.adjacent(x, 1);
    Box yAdj = domainBox.adjacent(y, 1);
    Box xEdge = domainBox.edge(x, 1);
    Box yEdge = domainBox.edge(y, 1);
    Box xyAdj = domainBox.adjacent(x+y,2);
    Box xyEdge = domainBox.edge(x+y,2);
    for (BlockIndex bi = 0; bi < domain.size(); bi++)
    {
        BlockIndex bx = (bi + 1) % XPOINT_SIZE;
        BlockIndex by = (bi + XPOINT_SIZE - 1) % XPOINT_SIZE;
        Point sx = x*domainSize*2;
        Point sy = y*domainSize*2;

        EXPECT_EQ(domain.convertNode(origin, bi, bx), sy);
        EXPECT_EQ(domain.convertNode(origin, bi, by), sx);
    
        EXPECT_EQ(domain.convertBox(xAdj, bi, bx), yEdge);
        EXPECT_EQ(domain.convertBox(yAdj, bi, by), xEdge);
        
        if (XPOINT_SIZE > 3)
        {
            BlockIndex bxy = (bi + 2) % XPOINT_SIZE;
            EXPECT_EQ(domain.convertBox(xyAdj, bi, bxy), xyEdge);
            EXPECT_EQ(domain.convertBox(xyEdge, bi, bxy), xyAdj);

        }
    }
}
#if DIM == 3
TEST(MBProblemDomain, CubedSphere)
{
    int domainSize = 4;
    int boxSize = 4;
    int thickness = 8;
    int radialDir = CUBED_SPHERE_SHELL_RADIAL_COORD;
    auto domain = CubedSphereShell::Domain(domainSize, thickness, radialDir);

}
#endif
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
