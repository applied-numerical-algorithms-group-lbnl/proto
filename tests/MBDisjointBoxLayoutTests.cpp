#include <gtest/gtest.h>
#include "Proto.H"
#include "TestFunctions.H"

using namespace Proto;

#if 1
TEST(MBDisjointBoxLayout, Iteration)
{
    int domainSize = 64;
    int boxSize = 16;
    int numBlocks = 5;
    auto domain = buildXPoint(domainSize, numBlocks);
    Point boxSizeVect = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);

    std::set<std::pair<unsigned int, Box>> correctPatchSet;
    Box patches = Box::Cube(domainSize).coarsen(boxSizeVect);
    for (unsigned int bi = 0; bi < domain.size(); bi++)
    {
        for (auto p : patches)
        {
            Box b(p, p);
            b = b.refine(boxSizeVect);
            std::pair<unsigned int, Box> patch(bi, b);
            correctPatchSet.insert(patch);
        }
    }

    int N = pow(domainSize / boxSize, DIM) * domain.size();
    int n0 = N / numProc();
    int r = N % numProc();
    if (procID() < r)
    {
        n0++;
    }
    std::set<std::pair<unsigned int, Box>> patchSet;
    for (auto index : layout)
    {
        Box b = layout[index];
        auto block = layout.block(index);
        std::pair<unsigned int, Box> patch(block, b);
        EXPECT_EQ(b.sizes(), boxSizeVect);
        EXPECT_EQ(patchSet.find(patch), patchSet.end());
        EXPECT_FALSE(correctPatchSet.find(patch) == correctPatchSet.end());
        for (auto pi : patchSet)
        {
            if (pi.first == block)
            {
                EXPECT_TRUE((pi.second & b).empty());
            }
        }
        patchSet.insert(patch);
    }
    EXPECT_EQ(patchSet.size(), n0);
}

TEST(MBDisjointBoxLayout, PatchConnectivity)
{
    int domainSize = 32;
    int boxSize = 16;
    int numBlocks = 5;
    auto domain = buildXPoint(domainSize, numBlocks);
    Point boxSizeVect = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);

    Box domainBox = Box::Cube(domainSize / boxSize);
    Box adjX = domainBox.adjacent(Point::X(), 1);
    Box adjY = domainBox.adjacent(Point::Y(), 1);
    Box adjXY = domainBox.adjacent(Point::X() + Point::Y(), 1);
    Box edgX = domainBox.edge(Point::X(), 1);
    Box edgY = domainBox.edge(Point::Y(), 1);
    Box edgXY = domainBox.edge(Point::X() + Point::Y(), 1);
    auto CW = CoordPermutation::cw();
    auto CCW = CoordPermutation::ccw();
    auto R = CW * CW;
    for (auto i1 : layout)
    {
        for (auto i2 : layout)
        {
            auto p1 = layout.patch(i1);
            auto p2 = layout.patch(i2);
            Point q2 = Point::Zeros();

            Point p = layout.connectivity(p1, p2);
            if (p1.block == p2.block)
            {
                q2 = p2.point;
            }
            else if (p2.block == (p1.block + 1) % layout.numBlocks())
            {
                q2 = CW.rotateCell(p2.point, edgY, adjX);
            }
            else if (p2.block == (p1.block + layout.numBlocks() - 1) % layout.numBlocks())
            {
                q2 = CCW.rotateCell(p2.point, edgX, adjY);
            }
            else
            {
                q2 = R.rotateCell(p2.point, edgXY, adjXY);
            }

            if (Box::Kernel(1).shift(p1.point).containsPoint(q2))
            {
                EXPECT_EQ(p, q2 - p1.point);
            }
            else
            {
                EXPECT_EQ(p, Point::Zeros());
            }
        }
    }
}

TEST(MBDisjointBoxLayout, Find)
{
    int domainSize = 64;
    int boxSize = 16;
    int numBlocks = 5;
    auto domain = buildXPoint(domainSize, numBlocks);
    Point boxSizeVect = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);

    for (auto iter : layout)
    {
        Point patch = layout.point(iter);
        auto block = layout.block(iter);
        auto index = layout.find(patch, block);
        EXPECT_EQ(index, iter);
    }
}

TEST(MBDisjointBoxLayout, Coarsen)
{
    int domainSize = 64;
    int boxSize = 16;
    int refRatio = 2;
    int numBlocks = 5;
    auto domain = buildXPoint(domainSize, numBlocks);
    Point boxSizeVect = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);

    std::vector<Point> refRatios(numBlocks, Point::Ones(refRatio));
    auto crseLayout = layout.coarsen(refRatios);

    auto crseDomain = buildXPoint(domainSize / refRatio);
    MBDisjointBoxLayout crseLayoutSoln(crseDomain, Point::Ones(boxSize / refRatio));

    EXPECT_EQ(crseLayout.numBoxes(), layout.numBoxes());
    EXPECT_TRUE(crseLayout.compatible(layout));
    EXPECT_TRUE(crseLayout.compatible(crseLayoutSoln));
    for (auto iter : crseLayout)
    {
        EXPECT_EQ(crseLayout[iter], crseLayoutSoln[iter]);
    }
}

TEST(MBDisjointBoxLayout, BoundaryQueries)
{
    for (int testNum = 0; testNum < 2; testNum++)
    {
        int domainSize = 16;
        int boxSize = 4;
        int numBlocks = 3 + 2 * testNum;
        auto domain = buildXPoint(domainSize, numBlocks);
        Point boxSizeVect = Point::Ones(boxSize);
        MBDisjointBoxLayout layout(domain, boxSizeVect);

        Box patchDomain = Box::Cube(domainSize / boxSize);
        Box blockBoundaryX = patchDomain.adjacent(Point::X(), 1);
        Box blockBoundaryY = patchDomain.adjacent(Point::Y(), 1);
        Box blockBoundaryXY = patchDomain.adjacent(Point::X() + Point::Y(), 1);
        for (auto iter : layout)
        {
            MBPoint patch = layout.patch(iter);
            BlockIndex block = layout.block(iter);
            BlockIndex adjBlockX = (block + 1) % numBlocks;
            BlockIndex adjBlockY = (block + numBlocks - 1) % numBlocks;
            for (auto dir : Point::Directions())
            {
                MBPoint adjPatch(patch.point + dir, patch.block);
                if (patchDomain.containsPoint(adjPatch.point))
                {
                    EXPECT_FALSE(layout.isBlockBoundary(patch, dir));
                    EXPECT_FALSE(layout.isDomainBoundary(patch, dir));
                    EXPECT_TRUE(layout.isInteriorBoundary(patch, dir));
                }
                else if (blockBoundaryX.containsPoint(adjPatch.point))
                {
                    EXPECT_TRUE(layout.isBlockBoundary(patch, dir));
                    for (auto bi = 0; bi < numBlocks; bi++)
                    {
                        EXPECT_EQ(layout.isBlockBoundary(patch, dir, bi), (bi == adjBlockX));
                    }
                    EXPECT_FALSE(layout.isDomainBoundary(patch, dir));
                    EXPECT_FALSE(layout.isInteriorBoundary(patch, dir));
                }
                else if (blockBoundaryY.containsPoint(adjPatch.point))
                {
                    EXPECT_TRUE(layout.isBlockBoundary(patch, dir));
                    for (auto bi = 0; bi < numBlocks; bi++)
                    {
                        EXPECT_EQ(layout.isBlockBoundary(patch, dir, bi), (bi == adjBlockY));
                    }
                    EXPECT_FALSE(layout.isDomainBoundary(patch, dir));
                    EXPECT_FALSE(layout.isInteriorBoundary(patch, dir));
                }
                else if (blockBoundaryXY.containsPoint(adjPatch.point))
                {
                    EXPECT_TRUE(layout.isBlockBoundary(patch, dir));
                    for (auto bi = adjBlockX + 1; bi % numBlocks != adjBlockY; bi++)
                    {
                        EXPECT_TRUE(layout.isBlockBoundary(patch, dir, bi % numBlocks));
                    }
                    EXPECT_FALSE(layout.isDomainBoundary(patch, dir));
                    EXPECT_FALSE(layout.isInteriorBoundary(patch, dir));
                }
                else
                {
                    EXPECT_FALSE(layout.isBlockBoundary(patch, dir));
                    EXPECT_TRUE(layout.isDomainBoundary(patch, dir));
                    EXPECT_FALSE(layout.isInteriorBoundary(patch, dir));
                }
            }
        }
    }
}

TEST(MBDisjointBoxLayout, PatchOnBoundaryQueries)
{
    int domainSize = 16;
    int boxSize = 4;
    int numBlocks = 5;
    auto domain = buildXPoint(domainSize, numBlocks);
    Point boxSizeVect = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);

    Box patchDomain = Box::Cube(domainSize / boxSize);
    auto X = Point::X();
    auto Y = Point::Y();
    for (auto iter : layout)
    {
        MBPoint patch = layout.patch(iter);
        BlockIndex block = layout.block(iter);
        BlockIndex adjBlockX = (block + 1) % numBlocks;
        BlockIndex adjBlockY = (block + numBlocks - 1) % numBlocks;
        if (patchDomain.grow(-1).containsPoint(patch.point))
        {
            for (BlockIndex bi = 0; bi < numBlocks; bi++)
            {
                EXPECT_FALSE(layout.isPatchOnBlockBoundary(patch));
                for (auto dir : Point::Directions())
                {
                    EXPECT_FALSE(layout.isPatchOnBlockBoundary(patch, dir));
                }
            }
        }
        else
        {
            for (auto dir : Point::Directions())
            {
                Box boundBox = patchDomain.edge(dir, 1);
                if (!boundBox.containsPoint(patch.point))
                {
                    continue;
                }
                if (dir == X)
                {
                    EXPECT_TRUE(layout.isPatchOnBlockBoundary(patch));
                    EXPECT_TRUE(layout.isPatchOnBlockBoundary(patch, X));
                }
                else if (dir == Y)
                {
                    EXPECT_TRUE(layout.isPatchOnBlockBoundary(patch));
                    EXPECT_TRUE(layout.isPatchOnBlockBoundary(patch, Y));
                }
                else if (dir == X + Y)
                {
                    EXPECT_TRUE(layout.isPatchOnBlockBoundary(patch));
                    EXPECT_TRUE(layout.isPatchOnBlockBoundary(patch, X + Y));
                }
                else
                {
                    EXPECT_TRUE(layout.isPatchOnDomainBoundary(patch));
                    EXPECT_TRUE(layout.isPatchOnDomainBoundary(patch, dir));
                }
            }
        }
    }
}

TEST(MBDisjointBoxLayout, PatchInBoundaryQueries)
{
    int domainSize = 16;
    int boxSize = 4;
    for (int testNum = 0; testNum < 3; testNum++)
    {
        int numBlocks = 3+testNum;
        auto domain = buildXPoint(domainSize, numBlocks);
        Point boxSizeVect = Point::Ones(boxSize);
        MBDisjointBoxLayout layout(domain, boxSizeVect);

        Box patchDomain = Box::Cube(domainSize / boxSize);
        auto X = Point::X();
        auto Y = Point::Y();

        for (BlockIndex block = 0; block < numBlocks; block++)
        {
            for (auto pi : patchDomain)
            {
                MBPoint patch(pi, block);
                EXPECT_FALSE(layout.isPatchInBlockBoundary(patch));
                EXPECT_FALSE(layout.isPatchInDomainBoundary(patch));
                EXPECT_FALSE(layout.isPatchInTriplePointRegion(patch));
            }
            for (auto dir : Point::Directions())
            {
                Box boundBox = patchDomain.adjacent(dir, 1);
                for (auto pi : boundBox)
                {
                    MBPoint patch(pi, block);
                    if (dir == X || dir == Y)
                    {
                        EXPECT_TRUE(layout.isPatchInBlockBoundary(patch));
                        EXPECT_TRUE(layout.isPatchInBlockBoundary(patch, dir));
                        EXPECT_FALSE(layout.isPatchInDomainBoundary(patch));
                        EXPECT_FALSE(layout.isPatchInTriplePointRegion(patch));
                    }
                    else if (dir == X + Y)
                    {
                        EXPECT_TRUE(layout.isPatchInBlockBoundary(patch));
                        EXPECT_TRUE(layout.isPatchInBlockBoundary(patch, dir));
                        EXPECT_FALSE(layout.isPatchInDomainBoundary(patch));
                        EXPECT_EQ(layout.isPatchInTriplePointRegion(patch), (numBlocks == 3));
                    }
                    else
                    {
                        EXPECT_FALSE(layout.isPatchInBlockBoundary(patch));
                        EXPECT_TRUE(layout.isPatchInDomainBoundary(patch));
                        EXPECT_TRUE(layout.isPatchInDomainBoundary(patch, dir));
                        EXPECT_FALSE(layout.isPatchInTriplePointRegion(patch));
                    }
                }
            }
        }
    }
}
#endif
TEST(MBDisjointBoxLayout, RefinementBoundaryQueries)
{
    for (int testNum = 0; testNum < 1; testNum++)
    {
        int domainSize = 16;
        int boxSize = 4;
        int numBlocks = 3 + 2 * testNum;

        auto domain = buildXPoint(domainSize, numBlocks);
        Point boxSizeVect = Point::Ones(boxSize);
        std::vector<MBPoint> patches;
        std::vector<Point> boxSizes;
        int patchesPerBlock1D = domainSize / boxSize - 1;
        Point skipPatch = (Point::X() + Point::Y())*patchesPerBlock1D;
        for (BlockIndex block = 0; block < numBlocks; block++)
        {
            boxSizes.push_back(boxSizeVect);
            for (auto patch : Box::Cube(domainSize / boxSize))
            {
                if (patch == skipPatch && block == 0) { continue; }
                patches.push_back(MBPoint(patch, block));
            }
        }
        MBDisjointBoxLayout layout(domain, patches, boxSizes);
        #if PR_VERBOSE > 0
        HDF5Handler h5;
        MBLevelBoxData<int, 1, HOST> data(layout, Point::Zeros());
        MBLevelMap<MBMap_XPointRigid<HOST>, HOST> map;
        map.define(layout, Point::Zeros());
        for (BlockIndex bi = 0; bi < numBlocks; bi++) { map[bi].setNumBlocks(numBlocks); }
        map.initialize();
        data.setVal(0);
        for (auto index : layout)
        {
            auto patch = layout.patch(index);
            if (layout.isPatchOnRefinementBoundary(patch))
            {
                auto& di = data[index];
                di.setVal(1);
            }
        }
        h5.writeMBLevel(map, data, "MBDBL_TestRefBounds_T%i", testNum);
        #endif
        for (auto index : layout)
        {
            auto patch = layout.patch(index);

            for (auto dir : Point::Directions())
            {
                Point adjPatch = patch.point + dir;
                if (layout.isDomainBoundary(patch, dir))
                {
                    EXPECT_FALSE(layout.isRefinementBoundary(patch, dir));
                } else if (layout.isBlockBoundary(patch, dir))
                {
                    adjPatch = layout.patchDomain().convertPoint(adjPatch, patch.block, 0, PR_CELL);
                    bool foundSkipPatch = (adjPatch == skipPatch);
                    foundSkipPatch &= (patch.block != 0);
                    bool isRefBound = layout.isRefinementBoundary(patch, dir);
                    EXPECT_EQ(isRefBound, foundSkipPatch);
                } else {
                    bool foundSkipPatch = (adjPatch == skipPatch);
                    foundSkipPatch &= (patch.block == 0);
                    bool isRefBound = layout.isRefinementBoundary(patch, dir);
                    EXPECT_EQ(isRefBound, foundSkipPatch);
                }
            }
        }
    }
}
#if 1
TEST(MBDisjointBoxLayout, Connectivity_XPoint)
{
    int domainSize = 8;
    int boxSize = 4;
    int numBlocks = 5;
    auto domain = buildXPoint(domainSize, numBlocks);
    Point boxSizeVect = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);

    Point patchDomainSize = Point::Ones(domainSize / boxSize);
    Box patchDomain(patchDomainSize);

    for (auto iter_i : layout)
    {
        for (auto iter_j : layout)
        {
            auto pi = layout.patch(iter_i);
            auto pj = layout.patch(iter_j);

            if (pi.block == 0 && pi.block != pj.block)
            {
                auto dij = layout.connectivity(pi, pj);
                auto dji = layout.connectivity(pj, pi);
                if (dij == Point::Zeros())
                {
                    EXPECT_EQ(dji, Point::Zeros());
                }
            }
        }
    }
}

#if DIM == 3
TEST(MBDisjointBoxLayout, Connectivity_CubedSphere)
{
    int domainSize = 8;
    int boxSize = 4;
    int thickness = 8;
    int radialDir = CUBED_SPHERE_SHELL_RADIAL_COORD;
    auto domain = CubedSphereShell::Domain(domainSize, thickness, radialDir);
    Point boxSizeVect = Point::Ones(boxSize);
    boxSizeVect[radialDir] = min(thickness, boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);

    Point patchDomainSize = Point::Ones(domainSize / boxSize);
    patchDomainSize[radialDir] = thickness > boxSize ? thickness / boxSize : 1;
    Box patchDomain(patchDomainSize);

    for (auto iter_i : layout)
    {
        for (auto iter_j : layout)
        {
            auto pi = layout.patch(iter_i);
            auto pj = layout.patch(iter_j);

            if (pi.block == 0 && pi.block != pj.block)
            {
                auto dij = layout.connectivity(pi, pj);
                auto dji = layout.connectivity(pj, pi);
                if (dij == Point::Zeros())
                {
                    EXPECT_EQ(dji, Point::Zeros());
                }
            }
        }
    }
}
#endif
#endif
int main(int argc, char *argv[])
{
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
