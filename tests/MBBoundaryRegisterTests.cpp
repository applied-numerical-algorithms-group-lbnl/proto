#include <gtest/gtest.h>
#include "Proto.H"
#include "Lambdas.H"

using namespace Proto;

TEST(MBBoundaryRegister, Construction) {
    int domainSize = 32;
    int boxSize = 16;
    int numBlocks = 5;
    int ghostSize = 0;
    int depth = 1;
    bool bothSides = true;
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

    for (int ti = 0; ti < 4; ti++)
    {
        switch (ti)
        {
            case 0: depth = +1; bothSides = true; break;
            case 1: depth = -1; bothSides = true; break;
            case 2: depth = +1; bothSides = false; break;
            case 3: depth = -1; bothSides = false; break;
        }
        MBBoundaryRegister<int, 1, HOST, PR_CELL> boundRegister(layout, depth, ghost, bothSides);
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
                auto bounds = boundRegister.bounds(iter, dir);
                if (bounds.size() > 0)
                {
                    //std::cout << "Checking boundaries from patch " << patchID << " in block " << blockID << " in direction " << dir << std::endl;
                }
                if (bothSides)
                {
                    for (auto bound : bounds)
                    {
                        EXPECT_EQ(
                                boundRegister.local(iter, bound.adjIndex).box(),
                                boundRegister.adjacent(iter, bound.adjIndex).box());
                    }
                }

                Box interior = patchBox.edge(dir, abs(depth));
                Box exterior = patchBox.adjacent(dir, abs(depth));
                Box intext = interior;
                intext &= exterior.low();
                intext &= exterior.high();

                Box patchBoundary;
                Box adjPatchBoundary;
                if (bothSides)
                {
                    patchBoundary = intext;
                    adjPatchBoundary = patchBoundary;
                }
                else if (depth < 0)
                {
                    patchBoundary = interior;
                    adjPatchBoundary = exterior;
                }
                else if (depth > 0)
                {
                    patchBoundary = exterior;
                    adjPatchBoundary = interior;
                }

                if (patchDomain.contains(neighbor))
                {
                    EXPECT_EQ(bounds.size(), 0);
                } else if (patchDomain.adjacent(nx,1).contains(neighbor))
                {
                    EXPECT_EQ(bounds.size(), 1);
                    EXPECT_TRUE(layout.isBlockBoundary(iter, dir, xBlock));
                    /*
                    Point adjDir = -CCW(dir);
                    Box adjPatchBoundary;
                    if (bothSides)
                    {
                        adjPatchBoundary = patchBoundary;
                        //adjPatchBoundary = adjPatchBox.edge(adjDir, abs(depth));
                        //adjPatchBoundary = adjPatchBoundary.extrude(adjDir, abs(depth));
                    } else if (depth < 0)
                    {
                        adjPatchBoundary = adjPatchBox.edge(adjDir, depth);
                    } else if (depth > 0)
                    {
                        adjPatchBoundary = adjPatchBox.adjacent(adjDir, depth);
                    }
                    */
                    EXPECT_EQ(layout.block(bounds[0].localIndex),  blockID);
                    EXPECT_EQ(layout.block(bounds[0].adjIndex),  xBlock);
                    EXPECT_EQ(bounds[0].localData->box(),  patchBoundary.grow(ghost));
                    EXPECT_EQ(bounds[0].adjData->box(),  adjPatchBoundary.grow(ghost));
                } else if (patchDomain.adjacent(ny,1).contains(neighbor))
                {
                    if (neighbor == Point::Y()*(domainSize/boxSize)){continue;} //we manually removed this patch. 
                    EXPECT_EQ(bounds.size(),  1);
                    EXPECT_TRUE(layout.isBlockBoundary(iter, dir, yBlock));
                    /*
                    Point adjDir = -CW(dir); 
                    Box adjPatchBoundary;
                    if (bothSides)
                    {
                        adjPatchBoundary = adjPatchBox.edge(adjDir, abs(depth));
                        adjPatchBoundary = adjPatchBoundary.extrude(adjDir, abs(depth));
                    } else if (depth < 0)
                    {
                        adjPatchBoundary = adjPatchBox.edge(adjDir, depth);
                    } else if (depth > 0)
                    {
                        adjPatchBoundary = adjPatchBox.adjacent(adjDir, depth);
                    }
                    */
                    EXPECT_EQ(layout.block(bounds[0].localIndex),  blockID);
                    EXPECT_EQ(layout.block(bounds[0].adjIndex),  yBlock);
                    EXPECT_EQ(bounds[0].localData->box(),  patchBoundary.grow(ghost));
                    EXPECT_EQ(bounds[0].adjData->box(),  adjPatchBoundary.grow(ghost));
                } else if (patchDomain.adjacent(nx+ny,1).contains(neighbor))
                {
                    EXPECT_EQ(bounds.size(), numBlocks-3);
                    /*
                    Point adjDir = -dir;
                    adjDir[0] = dir[0]; adjDir[1] = dir[1];
                    Box adjPatchBoundary;
                    if (bothSides)
                    {
                        adjPatchBoundary = adjPatchBox.edge(adjDir, abs(depth));
                        adjPatchBoundary = adjPatchBoundary.extrude(adjDir, abs(depth));
                    } else if (depth < 0)
                    {
                        adjPatchBoundary = adjPatchBox.edge(adjDir, depth);
                    } else if (depth > 0)
                    {
                        adjPatchBoundary = adjPatchBox.adjacent(adjDir, depth);
                    }
                    */
                    for (auto bound : bounds)
                    {
                        EXPECT_EQ(layout.block(bound.localIndex), blockID);
                        EXPECT_NE(layout.block(bound.adjIndex), blockID);
                        EXPECT_NE(layout.block(bound.adjIndex), yBlock);
                        EXPECT_NE(layout.block(bound.adjIndex), xBlock);
                        EXPECT_EQ(bound.localData->box(), patchBoundary.grow(ghost));
                        EXPECT_EQ(bound.adjData->box(), adjPatchBoundary.grow(ghost));
                    }
                } else {
                    EXPECT_EQ(bounds.size(), 0);
                }
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
