#include <gtest/gtest.h>
#include "Proto.H"
#include "Lambdas.H"

using namespace Proto;

TEST(MBBoundaryRegister, Construction) {
    int domainSize = 32;
    int boxSize = 16;
    int numBlocks = 5;
    int ghostSize = 1;
    int depth = 1;
    auto domain = buildXPoint(domainSize, numBlocks);
    std::vector<Point> boxSizeVect;
    std::vector<MBPatchID_t> patches;
    for (int bi = 0; bi < numBlocks; bi++)
    {
        boxSizeVect.push_back(Point::Ones(boxSize));
        for (auto pi : Box::Cube(domainSize / boxSize))
        {
            if (pi == Point::X()*(domainSize/boxSize - 1)) {continue;}
            patches.push_back(MBPatchID_t(pi,bi));
        }
    }
    MBDisjointBoxLayout layout(domain, patches, boxSizeVect);
    Point ghost = Point::Ones(ghostSize);
    MBBoundaryRegister<int, 1, HOST, PR_CELL> ghostRegister(layout, depth, ghost);
    MBBoundaryRegister<int, 1, HOST, PR_CELL> fluxRegister(layout, -depth, ghost);

    Point nx = Point::X();
    Point ny = Point::Y();
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
        unsigned int xBlock = (blockID+1) % numBlocks;
        unsigned int yBlock = (blockID-1+numBlocks) % numBlocks;
        auto blockLayout = layout.getBlock(blockID);
        Box patchBox = layout[iter]; 
        for (auto dir : K)
        {
            Point neighbor = patchID + dir;
            Point adjPatch = patchMap[neighbor];
            Box adjPatchBox = Box(adjPatch, adjPatch).refine(boxSize);
            auto bounds = ghostRegister.bounds(iter, dir);
            auto bounds2 = fluxRegister.bounds(iter, dir);
            
            if (patchDomain.contains(neighbor))
            {
                EXPECT_EQ(bounds.size(), 0);
                EXPECT_EQ(bounds2.size(), 0);
            } else if (patchDomain.adjacent(nx,1).contains(neighbor))
            {
                EXPECT_EQ(bounds.size(), 1);
                EXPECT_EQ(bounds2.size(), 1);
                EXPECT_TRUE(layout.isBlockBoundary(iter, dir, xBlock));
                Box patchBoundary  = patchBox.adjacent(dir, 1);
                Box patchBoundary2 = patchBox.edge(dir, 1);
                Point adjDir = -CCW(dir);
                Box adjPatchBoundary  = adjPatchBox.edge(adjDir, 1);
                Box adjPatchBoundary2 = adjPatchBox.adjacent(adjDir, 1);

                EXPECT_EQ(layout.block(bounds[0].localIndex),  blockID);
                EXPECT_EQ(layout.block(bounds2[0].localIndex), blockID);
                EXPECT_EQ(layout.block(bounds[0].adjIndex),  xBlock);
                EXPECT_EQ(layout.block(bounds2[0].adjIndex), xBlock);
                EXPECT_EQ(bounds[0].localData->box(),  patchBoundary.grow(ghost));
                EXPECT_EQ(bounds2[0].localData->box(), patchBoundary2.grow(ghost));
                EXPECT_EQ(bounds[0].adjData->box(),  adjPatchBoundary.grow(ghost));
                EXPECT_EQ(bounds2[0].adjData->box(), adjPatchBoundary2.grow(ghost));
            } else if (patchDomain.adjacent(ny,1).contains(neighbor))
            {
                if (neighbor == Point::Y()*(domainSize/boxSize)){continue;} //we manually removed this patch. 
                EXPECT_EQ(bounds.size(),  1);
                EXPECT_EQ(bounds2.size(), 1);
                EXPECT_TRUE(layout.isBlockBoundary(iter, dir, yBlock));
                Box patchBoundary   = patchBox.adjacent(dir, 1);
                Box patchBoundary2  = patchBox.edge(dir, 1);
                Point adjDir = -CW(dir); 
                Box adjPatchBoundary    = adjPatchBox.edge(adjDir, 1);
                Box adjPatchBoundary2   = adjPatchBox.adjacent(adjDir, 1);
                EXPECT_EQ(layout.block(bounds[0].localIndex),  blockID);
                EXPECT_EQ(layout.block(bounds2[0].localIndex), blockID);
                EXPECT_EQ(layout.block(bounds[0].adjIndex),  yBlock);
                EXPECT_EQ(layout.block(bounds2[0].adjIndex), yBlock);
                EXPECT_EQ(bounds[0].localData->box(),  patchBoundary.grow(ghost));
                EXPECT_EQ(bounds2[0].localData->box(), patchBoundary2.grow(ghost));
                EXPECT_EQ(bounds[0].adjData->box(),  adjPatchBoundary.grow(ghost));
                EXPECT_EQ(bounds2[0].adjData->box(), adjPatchBoundary2.grow(ghost));
            } else if (patchDomain.adjacent(nx+ny,1).contains(neighbor))
            {
                EXPECT_EQ(bounds.size(), numBlocks-3);
                EXPECT_EQ(bounds2.size(), numBlocks-3);
                Box patchBoundary = patchBox.adjacent(dir,1);
                Box patchBoundary2 = patchBox.edge(dir,1);
                Point adjDir = -dir;
                adjDir[0] = dir[0]; adjDir[1] = dir[1];
                Box adjPatchBoundary = adjPatchBox.edge(adjDir, 1);
                Box adjPatchBoundary2 = adjPatchBox.adjacent(adjDir, 1);
                for (auto bound : bounds)
                {
                    EXPECT_EQ(layout.block(bound.localIndex), blockID);
                    EXPECT_NE(layout.block(bound.adjIndex), blockID);
                    EXPECT_NE(layout.block(bound.adjIndex), yBlock);
                    EXPECT_NE(layout.block(bound.adjIndex), xBlock);
                    EXPECT_EQ(bound.localData->box(), patchBoundary.grow(ghost));
                    EXPECT_EQ(bound.adjData->box(), adjPatchBoundary.grow(ghost));
                }
                for (auto bound : bounds2)
                {
                    EXPECT_EQ(layout.block(bound.localIndex), blockID);
                    EXPECT_NE(layout.block(bound.adjIndex), blockID);
                    EXPECT_NE(layout.block(bound.adjIndex), yBlock);
                    EXPECT_NE(layout.block(bound.adjIndex), xBlock);
                    EXPECT_EQ(bound.localData->box(), patchBoundary2.grow(ghost));
                    EXPECT_EQ(bound.adjData->box(), adjPatchBoundary2.grow(ghost));
                }
            } else {
                EXPECT_EQ(bounds.size(), 0);
                EXPECT_EQ(bounds2.size(), 0);
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
