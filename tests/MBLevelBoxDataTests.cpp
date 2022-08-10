#include <gtest/gtest.h>
#include "ProtoMMB.H"
#include "Lambdas.H"

#define NCOMP 1
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
    auto domain = buildXPoint(domainSize);
    Point boxSizeVect = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);
    std::array<Point, DIM+1> ghost;
    ghost.fill(Point::Ones());
    MBLevelBoxData<int, NCOMP, HOST> hostData(layout, ghost);

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
            //Point ghostDir = dir*ghost;
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
                EXPECT_TRUE(hostData.isBlockBoundary(iter, dir, xBlock));
                Box patchBoundary = patchBox.adjacent(dir, 1);
                Point adjDir = -CCW(dir);
                //Point adjGhostDir = -CCW(ghostDir); 
                //Box adjPatchBoundary = adjPatchBox.edge(adjGhostDir);
                Box adjPatchBoundary = adjPatchBox.edge(adjDir, 1);

                EXPECT_EQ(layout.block(bounds[0].localIndex), blockID);
                EXPECT_EQ(layout.block(bounds[0].adjIndex), xBlock);
                EXPECT_EQ(bounds[0].localData->box(), patchBoundary);
                EXPECT_EQ(bounds[0].adjData->box(), adjPatchBoundary);
            } else if (patchDomain.adjacent(ny,1).contains(neighbor))
            {
                EXPECT_EQ(bounds.size(), 1);
                EXPECT_TRUE(hostData.isBlockBoundary(iter, dir, yBlock));
                Box patchBoundary = patchBox.adjacent(dir, 1);
                //Point adjGhostDir = -CW(ghostDir); 
                //Box adjPatchBoundary = adjPatchBox.edge(adjGhostDir);
                Point adjDir = -CW(dir); 
                Box adjPatchBoundary = adjPatchBox.edge(adjDir, 1);
                EXPECT_EQ(layout.block(bounds[0].localIndex), blockID);
                EXPECT_EQ(layout.block(bounds[0].adjIndex), yBlock);
                EXPECT_EQ(bounds[0].localData->box(), patchBoundary);
                EXPECT_EQ(bounds[0].adjData->box(), adjPatchBoundary);
            } else if (patchDomain.adjacent(nx+ny,1).contains(neighbor))
            {
                EXPECT_EQ(bounds.size(), XPOINT_SIZE-3);
                Box patchBoundary = patchBox.adjacent(dir,1);
                //Point adjGhostDir = -ghostDir;
                //adjGhostDir[0] = ghostDir[0]; adjGhostDir[1] = ghostDir[1];
                Point adjDir = -dir;
                adjDir[0] = dir[0]; adjDir[1] = dir[1];
                Box adjPatchBoundary = adjPatchBox.edge(adjDir, 1);
                for (auto bound : bounds)
                {
                    EXPECT_EQ(layout.block(bound.localIndex), blockID);
                    EXPECT_NE(layout.block(bound.adjIndex), blockID);
                    EXPECT_NE(layout.block(bound.adjIndex), yBlock);
                    EXPECT_NE(layout.block(bound.adjIndex), xBlock);
                    EXPECT_EQ(bound.localData->box(), patchBoundary);
                    EXPECT_EQ(bound.adjData->box(), adjPatchBoundary);
                }
            } else {
                EXPECT_EQ(bounds.size(), 0);
            }
        }
    }
}

TEST(MBLevelBoxData, Initialization) {
    int domainSize = 64;
    int boxSize = 16;
    int ghostSize = 1;
    auto domain = buildXPoint(domainSize);
    Point boxSizeVect = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);

    MBLevelBoxData<int, NCOMP, HOST> hostData(layout, Point::Ones(ghostSize));
    hostData.initialize(f_MBPointID);

    for (auto iter : layout)
    {
        Box patchBox = layout[iter];
        auto& hostData_i = hostData[iter];
        int block = layout.block(iter);
        BoxData<int, NCOMP, HOST> tempData(patchBox);
        forallInPlace_p(f_MBPointID, tempData, block);
        for (auto pi : patchBox)
        {
            for (int cc = 0; cc < NCOMP; cc++)
            {
                EXPECT_EQ(tempData(pi, cc), hostData_i(pi, cc));
            }
        }
    }
}

TEST(MBLevelBoxData, FillBoundaries) {
    int domainSize = 32;
    int boxSize = 16;
    int ghostSize = 1;
    auto domain = buildXPoint(domainSize);
    Point boxSizeVect = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);

    MBLevelBoxData<int, NCOMP, HOST> hostData(layout, Point::Ones(ghostSize));
    hostData.initialize(f_MBPointID);
    hostData.fillBoundaries();

    Box dirs = Box::Kernel(1);

    for (auto iter : layout)
    {
        auto block = layout.block(iter);
        auto patch = layout.point(iter);
        for (auto dir : dirs)
        {
            std::cout << "block: " << block << " | patch: " << patch << " | dir: " << dir << std::endl;
            auto bounds = hostData.bounds(iter, dir);
            for (auto bound : bounds)
            {
                std::cout << "localData: " << std::endl;
                bound.localData->printData();
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
