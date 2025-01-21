#include <gtest/gtest.h>
#include "Proto.H"
#include "TestFunctions.H"


using namespace Proto;

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

#if 1
TEST(MBProblemDomain, BoundaryQueries)
{
    int domainSize = 4;
    int ghostSize = 2;

    auto X = Point::X();
    auto Y = Point::Y();

    std::set<Point> blockBoundaries;
    blockBoundaries.insert(X);
    blockBoundaries.insert(Y);
    blockBoundaries.insert(X + Y);
    std::set<Point> domainBoundaries;
    for (auto dir : Box::Kernel(1))
    {
        if (dir == Point::Zeros())
        {
            continue;
        }
        if (blockBoundaries.count(dir) > 0)
        {
            continue;
        }
        domainBoundaries.insert(dir);
    }

    for (int testNum = 0; testNum < 2; testNum++)
    {
        int numBlocks;
        switch (testNum)
        {
            case 0: numBlocks = 5; break;
            case 1: numBlocks = 3; break;
        }
        auto domain = buildXPoint(domainSize, numBlocks);

        std::set<Point> triplePointBoundaries;
        if (numBlocks == 3)
        {
            triplePointBoundaries.insert(X + Y);
        }

        Box domainBox = Box::Cube(domainSize);
        Box searchBox = domainBox.grow(ghostSize);

        for (BlockIndex block = 0; block < numBlocks; block++)
        {
            for (auto point : searchBox)
            {
                // std::cout << "\n==================================================================" << std::endl;
                // std::cout << "checking point: " << point << " in block: " << block << std::endl;
                if (domainBox.containsPoint(point))
                {
                    EXPECT_TRUE(domain.isPointInInterior(point, block));
                    EXPECT_FALSE(domain.isPointInBlockBoundary(point, block));
                    EXPECT_FALSE(domain.isPointInDomainBoundary(point, block));
                    EXPECT_FALSE(domain.isPointInTriplePoint(point, block));
                }
                for (auto dir : Box::Kernel(1))
                {
                    if (dir == Point::Zeros()) { continue; }
                    if (domainBox.adjacent(dir, ghostSize).containsPoint(point))
                    {
                        // std::cout << "\tPoint is in boundary region with dir = " << dir << std::endl;
                        if (blockBoundaries.count(dir) > 0)
                        {
                            EXPECT_FALSE(domain.isPointInInterior(point, block));
                            EXPECT_TRUE(domain.isPointInBlockBoundary(point, block, dir));
                            EXPECT_FALSE(domain.isPointInDomainBoundary(point, block, dir));
                        }
                        else
                        {
                            EXPECT_FALSE(domain.isPointInBlockBoundary(point, block, dir));
                        }
                        if (domainBoundaries.count(dir) > 0)
                        {
                            EXPECT_FALSE(domain.isPointInInterior(point, block));
                            EXPECT_FALSE(domain.isPointInBlockBoundary(point, block, dir));
                            EXPECT_TRUE(domain.isPointInDomainBoundary(point, block, dir));
                        }
                        else
                        {
                            EXPECT_FALSE(domain.isPointInDomainBoundary(point, block, dir));
                        }
                        if (triplePointBoundaries.count(dir) > 0)
                        {
                            EXPECT_TRUE(domain.isPointInBlockBoundary(point, block, dir));
                            EXPECT_TRUE(domain.isPointInTriplePoint(point, block, dir));
                            EXPECT_FALSE(domain.isPointInDomainBoundary(point, block, dir));
                        }
                        else
                        {
                            EXPECT_FALSE(domain.isPointInTriplePoint(point, block, dir));
                        }
                    } else {
                        EXPECT_FALSE(domain.isPointInBlockBoundary(point, block, dir));
                        EXPECT_FALSE(domain.isPointInTriplePoint(point, block, dir));
                        EXPECT_FALSE(domain.isPointInDomainBoundary(point, block, dir));
                    }
                }
            }
        }

    }
}
#endif
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
    int numBlocks = 5;
    auto domain = buildXPoint(domainSize, numBlocks);
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
        BlockIndex bx = (bi + 1) % numBlocks;
        BlockIndex by = (bi + numBlocks - 1) % numBlocks;
        Point sx = x*domainSize*2;
        Point sy = y*domainSize*2;

        EXPECT_EQ(domain.convertNode(origin, bi, bx), sy);
        EXPECT_EQ(domain.convertNode(origin, bi, by), sx);
    
        EXPECT_EQ(domain.convertBox(xAdj, bi, bx), yEdge);
        EXPECT_EQ(domain.convertBox(yAdj, bi, by), xEdge);
        
        if (numBlocks > 3)
        {
            BlockIndex bxy = (bi + 2) % numBlocks;
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
