#include <gtest/gtest.h>
#include "Proto.H"
#include "TestFunctions.H"

using namespace Proto;
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
            auto b1 = layout.block(i1);
            auto b2 = layout.block(i2);
            auto p1 = layout.point(i1);
            auto p2 = layout.point(i2);
            Point q2 = Point::Zeros();

            Point p = layout.connectivity(i1, i2);
            if (b1 == b2)
            {
                q2 = p2;
            }
            else if (b2 == (b1 + 1) % layout.numBlocks())
            {
                q2 = CW.rotateCell(p2, edgY, adjX);
            }
            else if (b2 == (b1 + layout.numBlocks() - 1) % layout.numBlocks())
            {
                q2 = CCW.rotateCell(p2, edgX, adjY);
            }
            else
            {
                q2 = R.rotateCell(p2, edgXY, adjXY);
            }

            if (Box::Kernel(1).shift(p1).containsPoint(q2))
            {
                EXPECT_EQ(p, q2 - p1);
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
            PatchID patch = layout.point(iter);
            BlockIndex block = layout.block(iter);
            BlockIndex adjBlockX = (block + 1) % numBlocks;
            BlockIndex adjBlockY = (block + numBlocks - 1) % numBlocks;
            for (auto dir : Point::Directions())
            {
                PatchID adjPatch = patch + dir;
                if (patchDomain.containsPoint(adjPatch))
                {
                    EXPECT_FALSE(layout.isBlockBoundary(iter, dir));
                    EXPECT_FALSE(layout.isDomainBoundary(iter, dir));
                    EXPECT_TRUE(layout.isInteriorBoundary(iter, dir));
                }
                else if (blockBoundaryX.containsPoint(adjPatch))
                {
                    EXPECT_TRUE(layout.isBlockBoundary(iter, dir));
                    for (auto bi = 0; bi < numBlocks; bi++)
                    {
                        EXPECT_EQ(layout.isBlockBoundary(iter, dir, bi), (bi == adjBlockX));
                    }
                    EXPECT_FALSE(layout.isDomainBoundary(iter, dir));
                    EXPECT_FALSE(layout.isInteriorBoundary(iter, dir));
                }
                else if (blockBoundaryY.containsPoint(adjPatch))
                {
                    EXPECT_TRUE(layout.isBlockBoundary(iter, dir));
                    for (auto bi = 0; bi < numBlocks; bi++)
                    {
                        EXPECT_EQ(layout.isBlockBoundary(iter, dir, bi), (bi == adjBlockY));
                    }
                    EXPECT_FALSE(layout.isDomainBoundary(iter, dir));
                    EXPECT_FALSE(layout.isInteriorBoundary(iter, dir));
                }
                else if (blockBoundaryXY.containsPoint(adjPatch))
                {
                    EXPECT_TRUE(layout.isBlockBoundary(iter, dir));
                    for (auto bi = adjBlockX + 1; bi % numBlocks != adjBlockY; bi++)
                    {
                        EXPECT_TRUE(layout.isBlockBoundary(iter, dir, bi % numBlocks));
                    }
                    EXPECT_FALSE(layout.isDomainBoundary(iter, dir));
                    EXPECT_FALSE(layout.isInteriorBoundary(iter, dir));
                }
                else
                {
                    EXPECT_FALSE(layout.isBlockBoundary(iter, dir));
                    EXPECT_TRUE(layout.isDomainBoundary(iter, dir));
                    EXPECT_FALSE(layout.isInteriorBoundary(iter, dir));
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
        PatchID patch = layout.point(iter);
        BlockIndex block = layout.block(iter);
        BlockIndex adjBlockX = (block + 1) % numBlocks;
        BlockIndex adjBlockY = (block + numBlocks - 1) % numBlocks;
        if (patchDomain.grow(-1).containsPoint(patch))
        {
            for (BlockIndex bi = 0; bi < numBlocks; bi++)
            {
                EXPECT_FALSE(layout.isPatchOnBlockBoundary(patch, block));
                for (auto dir : Point::Directions())
                {
                    EXPECT_FALSE(layout.isPatchOnBlockBoundary(patch, block, dir));
                }
            }
        }
        else
        {
            for (auto dir : Point::Directions())
            {
                Box boundBox = patchDomain.edge(dir, 1);
                if (!boundBox.containsPoint(patch))
                {
                    continue;
                }
                if (dir == X)
                {
                    EXPECT_TRUE(layout.isPatchOnBlockBoundary(patch, block));
                    EXPECT_TRUE(layout.isPatchOnBlockBoundary(patch, block, X));
                }
                else if (dir == Y)
                {
                    EXPECT_TRUE(layout.isPatchOnBlockBoundary(patch, block));
                    EXPECT_TRUE(layout.isPatchOnBlockBoundary(patch, block, Y));
                }
                else if (dir == X + Y)
                {
                    EXPECT_TRUE(layout.isPatchOnBlockBoundary(patch, block));
                    EXPECT_TRUE(layout.isPatchOnBlockBoundary(patch, block, X + Y));
                }
                else
                {
                    EXPECT_TRUE(layout.isPatchOnDomainBoundary(patch, block));
                    EXPECT_TRUE(layout.isPatchOnDomainBoundary(patch, block, dir));
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
            for (auto patch : patchDomain)
            {
                EXPECT_FALSE(layout.isPatchInBlockBoundary(patch, block));
                EXPECT_FALSE(layout.isPatchInDomainBoundary(patch, block));
                EXPECT_FALSE(layout.isPatchInTriplePointRegion(patch, block));
            }
            for (auto dir : Point::Directions())
            {
                Box boundBox = patchDomain.adjacent(dir, 1);
                for (auto patch : boundBox)
                {
                    if (dir == X || dir == Y)
                    {
                        EXPECT_TRUE(layout.isPatchInBlockBoundary(patch, block));
                        EXPECT_TRUE(layout.isPatchInBlockBoundary(patch, block, dir));
                        EXPECT_FALSE(layout.isPatchInDomainBoundary(patch, block));
                        EXPECT_FALSE(layout.isPatchInTriplePointRegion(patch, block));
                    }
                    else if (dir == X + Y)
                    {
                        EXPECT_TRUE(layout.isPatchInBlockBoundary(patch, block));
                        EXPECT_TRUE(layout.isPatchInBlockBoundary(patch, block, dir));
                        EXPECT_FALSE(layout.isPatchInDomainBoundary(patch, block));
                        EXPECT_EQ(layout.isPatchInTriplePointRegion(patch, block), (numBlocks == 3));
                    }
                    else
                    {
                        EXPECT_FALSE(layout.isPatchInBlockBoundary(patch, block));
                        EXPECT_TRUE(layout.isPatchInDomainBoundary(patch, block));
                        EXPECT_TRUE(layout.isPatchInDomainBoundary(patch, block, dir));
                        EXPECT_FALSE(layout.isPatchInTriplePointRegion(patch, block));
                    }
                }
            }
        }
    }
}

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
            auto pi = layout.point(iter_i);
            auto bi = layout.block(iter_i);
            auto pj = layout.point(iter_j);
            auto bj = layout.block(iter_j);

            if (bi == 0 && bi != bj)
            {
                auto dij = layout.connectivity(iter_i, iter_j);
                auto dji = layout.connectivity(iter_j, iter_i);
                if (dij == Point::Zeros())
                {
                    EXPECT_EQ(dji, Point::Zeros());
                }
            }
        }
    }
}
#if 1
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
            auto pi = layout.point(iter_i);
            auto bi = layout.block(iter_i);
            auto pj = layout.point(iter_j);
            auto bj = layout.block(iter_j);

            if (bi == 0 && bi != bj)
            {
                auto dij = layout.connectivity(iter_i, iter_j);
                auto dji = layout.connectivity(iter_j, iter_i);
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
