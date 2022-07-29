#include <gtest/gtest.h>
#include "ProtoMMB.H"

#define NCOMP 3
#define XPOINT_SIZE 5
using namespace Proto;

MBProblemDomain buildXPoint(int a_domainSize)
{
    MBProblemDomain domain(XPOINT_SIZE);
    auto CCW = CoordPermutation::ccw();
    for (int ii = 0; ii < XPOINT_SIZE; ii++)
    {
        domain.defineBoundary(ii, (ii+1) % XPOINT_SIZE, 0, Side::Hi, CCW);
    }
    for (int bi = 0; bi < XPOINT_SIZE; bi++)
    {
        domain.defineDomain(bi, Point::Ones(a_domainSize));
    }
    return domain;
}

TEST(MBLevelBoxData, Construction) {
    int domainSize = 64;
    int boxSize = 16;
    int ghostSize = 1;
    auto domain = buildXPoint(domainSize);
    Point boxSizeVect = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);

    MBLevelBoxData<int, NCOMP, HOST> hostData(layout, Point::Ones(ghostSize));
   
    Point nx = Point::Basis(0);
    Point ny = Point::Basis(1);
    Box K = Box::Kernel(1);
    Box patchDomain = Box::Cube(domainSize).coarsen(boxSize);
    
    std::map<Point, Point> patchMap;
    for (auto pi : patchDomain.adjacent(nx,1))
    {
        Point p = pi;
        p[0] = pi[1];
        p[1] = patchDomain.size(1)-1;
        patchMap[pi] = p;
    }
    for (auto pi : patchDomain.adjacent(ny,1))
    {
        Point p = pi;
        p[1] = pi[0];
        p[0] = patchDomain.size(0)-1;
        patchMap[pi] = p;
    }
    for (auto pi : patchDomain.adjacent(nx+ny,1))
    {
        patchMap[pi] = pi-(nx+ny);
    }

    auto CW = CoordPermutation::cw();
    auto CCW = CoordPermutation::ccw();

    for (auto iter : layout)
    {
        auto patchID = layout.point(iter);
        auto blockID = layout.block(iter);

        unsigned int xBlock = (blockID+1) % XPOINT_SIZE;
        unsigned int yBlock = (blockID-1+XPOINT_SIZE) % XPOINT_SIZE;
        auto blockLayout = layout.blockLayout(blockID);
        Box patchBox = layout[iter]; 
        for (auto dir : K)
        {
            Point neighbor = patchID + dir;
            Point adjPatch = patchMap[neighbor];
            Box adjPatchBox = Box(adjPatch, adjPatch).refine(boxSize);
            auto bounds = hostData.bounds(iter, dir);
            if (patchDomain.contains(neighbor))
            {
                EXPECT_EQ(bounds.size(), 0);
            } else if (patchDomain.adjacent(nx,1).contains(neighbor))
            {
                EXPECT_EQ(bounds.size(), 1);
                Box patchBoundary = patchBox.adjacent(dir,1);
                Point adjDir = -CCW(dir); 
                Box adjPatchBoundary = adjPatchBox.edge(adjDir,1);

                EXPECT_EQ(bounds[0].localBlock, blockID);
                EXPECT_EQ(bounds[0].adjBlock, xBlock);
                EXPECT_EQ(bounds[0].localData->box(), patchBoundary);
                EXPECT_EQ(bounds[0].adjData->box(), adjPatchBoundary);
            } else if (patchDomain.adjacent(ny,1).contains(neighbor))
            {
                EXPECT_EQ(bounds.size(), 1);
                Box patchBoundary = patchBox.adjacent(dir,1);
                Point adjDir = -CW(dir); 
                Box adjPatchBoundary = adjPatchBox.edge(adjDir,1);
                EXPECT_EQ(bounds[0].localBlock, blockID);
                EXPECT_EQ(bounds[0].adjBlock, yBlock);
                EXPECT_EQ(bounds[0].localData->box(), patchBoundary);
                EXPECT_EQ(bounds[0].adjData->box(), adjPatchBoundary);
            } else if (patchDomain.adjacent(nx+ny,1).contains(neighbor))
            {
                EXPECT_EQ(bounds.size(), XPOINT_SIZE-3);
                Box patchBoundary = patchBox.adjacent(dir,1);
                Point adjDir = -dir;
                adjDir[0] = dir[0]; adjDir[1] = dir[1];
                Box adjPatchBoundary = adjPatchBox.edge(adjDir,1);
                for (auto bound : bounds)
                {
                    EXPECT_EQ(bound.localBlock, blockID);
                    EXPECT_NE(bound.adjBlock, blockID);
                    EXPECT_NE(bound.adjBlock, yBlock);
                    EXPECT_NE(bound.adjBlock, xBlock);
                    EXPECT_EQ(bound.localData->box(), patchBoundary);
                    EXPECT_EQ(bound.adjData->box(), adjPatchBoundary);
                }
            } else {
                EXPECT_EQ(bounds.size(), 0);
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
