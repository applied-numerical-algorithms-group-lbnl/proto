#include <gtest/gtest.h>
#include "Proto.H"
#include "TestFunctions.H"
using namespace Proto;

namespace {
    MBDisjointBoxLayout testLayout(int boxSize, int domainSize, int numBlocks)
    {
        auto domain = buildXPoint(domainSize, numBlocks);
        // std::vector<MBPoint> patches;
        // std::vector<Point> boxSizes;
        // for (BlockIndex bi = 0; bi < numBlocks; bi++)
        // {
        //     for (auto pi : Box::Cube(domainSize / boxSize))
        //     {
        //         if (pi != Point::X())
        //         {
        //             patches.push_back(MBPoint(pi, bi));
        //         }
        //     }
        //     boxSizes.push_back(Point::Ones(boxSize));
        // }
        //return MBDisjointBoxLayout(domain, patches, boxSizes);
        return MBDisjointBoxLayout(domain, Point::Ones(boxSize));
    }

    BoxData<int> printInterpLayout(MBInterpLayout layout, Point center, Point boundaryDir, MBIndex index)
    {
        auto footprint = layout.footprint(center, boundaryDir, index);
        Box b;
        for (auto fi : footprint)
        {
            b &= (fi + center);
        }
        BoxData<int> data(b);
        data.setVal(0);
        for (auto fi : footprint)
        {
            data(fi + center) = 1.0;
        }
        data(center) = 2.0;
        return data;
    }
}
TEST(MBInterpLayout, BaseFootprint)
{
    int boxSize = 8;
    int domainSize = boxSize;
    int ghostSize = 3;
    int numBlocks = 3;

    auto layout = testLayout(boxSize, domainSize, numBlocks);
    MBInterpLayout interpLayout(layout);
    std::set<Point> baseFootprint;
    for (auto pi : Box::Kernel(2))
    {
        if (pi.sumAbs() <= 2)
        {
            interpLayout.addPoint(pi);
            baseFootprint.insert(pi);
        }
    }
    for (auto pi : interpLayout.baseFootprint())
    {
        EXPECT_TRUE(baseFootprint.count(pi) == 1);
    }
}
TEST(MBInterpLayout, BlockBoundary)
{
    int boxSize = 8;
    int domainSize = boxSize;
    int ghostSize = 3;
    int numBlocks = 3;

    auto layout = testLayout(boxSize, domainSize, numBlocks);
    MBInterpLayout interpLayout(layout);
    std::set<Point> baseFootprint;
    for (auto pi : Box::Kernel(2))
    {
        if (pi.sumAbs() <= 2)
        {
            interpLayout.addPoint(pi);
            baseFootprint.insert(pi);
        }
    }
    Box domainBox = Point::Ones(domainSize);
    auto X = Point::X(); auto Y = Point::Y();
    Box yBound = domainBox.adjacent(Y, ghostSize);
    yBound = yBound.grow((Point::Ones() - Y)*(-2));
    Box xBound = domainBox.adjacent(X, ghostSize);
    xBound = yBound.grow((Point::Ones() - X)*(-2));

    for (auto iter : layout)
    {
        auto localXBound = layout[iter].adjacent(X, ghostSize);
        for (auto center : xBound & localXBound)
        {
            auto footprint = interpLayout.footprint(center, X, iter);
            for (auto fi : baseFootprint)
            {
                EXPECT_TRUE(footprint.count(fi) == 1);
            }
        }
        auto localYBound = layout[iter].adjacent(Y, ghostSize);
        for (auto center : yBound & localYBound)
        {
            auto footprint = interpLayout.footprint(center, Y, iter);
            for (auto fi : baseFootprint)
            {
                EXPECT_TRUE(footprint.count(fi) == 1);
            }
        }
    }
}
TEST(MBInterpLayout, DomainBoundary)
{
    int boxSize = 8;
    int domainSize = boxSize;
    int ghostSize = 3;
    int numBlocks = 3;

    auto layout = testLayout(boxSize, domainSize, numBlocks);
    MBInterpLayout interpLayout(layout);
    std::set<Point> baseFootprint;
    for (auto pi : Box::Kernel(2))
    {
        if (pi.sumAbs() <= 2)
        {
            interpLayout.addPoint(pi);
            baseFootprint.insert(pi);
        }
    }
    auto X = Point::X(); auto Y = Point::Y();
    std::set<Point> yBoundFootprint0;
    std::set<Point> yBoundFootprint1;
    std::set<Point> xBoundFootprint0;
    std::set<Point> xBoundFootprint1;
    for (auto pi : baseFootprint)
    {
        if (pi[1] >= 0)
        {
            xBoundFootprint0.insert(pi);
            xBoundFootprint0.insert(pi + Y);
        }
        if (pi[0] >= 0)
        {
            yBoundFootprint0.insert(pi);
            yBoundFootprint0.insert(pi + X);
        }
        if (pi != X*(-2)) { yBoundFootprint1.insert(pi); }
        if (pi != Y*(-2)) { xBoundFootprint1.insert(pi); }
    }

    Box domainBox = Point::Ones(domainSize);
    Box yBound0 = domainBox.adjacent(Y, ghostSize).edge(-X);
    Box yBound1 = yBound0.adjacent(X);
    Box xBound0 = domainBox.adjacent(X, ghostSize).edge(-Y);
    Box xBound1 = xBound0.adjacent(Y);
    for (auto dd = 2; dd < DIM; dd++)
    {
        yBound0 = yBound0.grow(dd,-2);
        yBound1 = yBound1.grow(dd,-2);
        xBound0 = xBound0.grow(dd,-2);
        xBound1 = xBound1.grow(dd,-2);
    }
    
    for (auto iter : layout)
    {
        auto localXBound = layout[iter].adjacent(X, ghostSize);
        for (auto center : xBound0 & localXBound)
        {
            auto footprint = interpLayout.footprint(center, X, iter);
            for (auto fi : xBoundFootprint0)
            {
                EXPECT_TRUE(footprint.count(fi) == 1);
            }
        }
        for (auto center : xBound1 & localXBound)
        {
            auto footprint = interpLayout.footprint(center, X, iter);
            for (auto fi : xBoundFootprint1)
            {
                EXPECT_TRUE(footprint.count(fi) == 1);
            }
        }
        auto localYBound = layout[iter].adjacent(Y, ghostSize);
        for (auto center : yBound0 & localYBound)
        {
            auto footprint = interpLayout.footprint(center, Y, iter);
            for (auto fi : yBoundFootprint0)
            {
                EXPECT_TRUE(footprint.count(fi) == 1);
            }
        }
        for (auto center : yBound1 & localYBound)
        {
            auto footprint = interpLayout.footprint(center, Y, iter);
            for (auto fi : yBoundFootprint1)
            {
                EXPECT_TRUE(footprint.count(fi) == 1);
            }
        }
    }
}
TEST(MBInterpLayout, TriplePointAdjacent)
{
    int boxSize = 8;
    int domainSize = boxSize;
    int ghostSize = 3;
    int numBlocks = 3;

    auto layout = testLayout(boxSize, domainSize, numBlocks);
    MBInterpLayout interpLayout(layout);
    std::set<Point> baseFootprint;
    for (auto pi : Box::Kernel(2))
    {
        if (pi.sumAbs() <= 2)
        {
            interpLayout.addPoint(pi);
            baseFootprint.insert(pi);
        }
    }
    
    auto X = Point::X(); auto Y = Point::Y();
    Box domainBox = Box::Cube(domainSize);
    Box yBoundBox = domainBox.adjacent(Y, ghostSize).edge(X);
    Box xBoundBox = domainBox.adjacent(X, ghostSize).edge(Y);
    for (int dd = 2; dd < DIM; dd++)
    {
        yBoundBox = yBoundBox.grow(dd,-2);
        xBoundBox = xBoundBox.grow(dd,-2);
    }

    for (auto iter : layout)
    {
        Box localXBoundBox = layout[iter].adjacent(X, ghostSize);
        for (auto center : localXBoundBox & xBoundBox)
        {
            auto footprint = interpLayout.footprint(center, X, iter);
            std::set<Point> footprintSln;
            for (auto fi : baseFootprint)
            {
                if (fi[1] <= 0)
                {
                    footprintSln.insert(fi);
                    Point corner = (X + Y)*(domainSize - 1);
                    Point p0 = center + fi - corner;
                    Point p1 = p0;
                    p1[0] = p0[1]; p1[1] = p0[0];
                    p1 += corner;
                    footprintSln.insert(p1 - center);
                }
            }

            for (auto fi : footprint)
            {
                EXPECT_TRUE(footprintSln.count(fi) == 1);
            }
        }
        Box localYBoundBox = layout[iter].adjacent(Y, ghostSize);
        for (auto center : localYBoundBox & yBoundBox)
        {
            auto footprint = interpLayout.footprint(center, Y, iter);
            std::set<Point> footprintSln;
            for (auto fi : baseFootprint)
            {
                if (fi[0] <= 0)
                {
                    footprintSln.insert(fi);
                    Point corner = (X + Y)*(domainSize - 1);
                    Point p0 = center + fi - corner;
                    Point p1 = p0;
                    p1[0] = p0[1]; p1[1] = p0[0];
                    p1 += corner;
                    footprintSln.insert(p1 - center);
                }
            }

            for (auto fi : footprint)
            {
                EXPECT_TRUE(footprintSln.count(fi) == 1);
            }
        }
    }

}

TEST(MBInterpLayout, TriplePointRegion)
{
    int boxSize = 8;
    int domainSize = boxSize;
    int ghostSize = 3;
    int numBlocks = 3;

    auto layout = testLayout(boxSize, domainSize, numBlocks);
    MBInterpLayout interpLayout(layout);
    std::set<Point> baseFootprint;
    for (auto pi : Box::Kernel(2))
    {
        if (pi.sumAbs() <= 2)
        {
            interpLayout.addPoint(pi);
            baseFootprint.insert(pi);
        }
    }
    
    auto X = Point::X(); auto Y = Point::Y();
    Box domainBox = Box::Cube(domainSize);
    Box triplePointRegion = domainBox.adjacent(X + Y, ghostSize);
    Box shrunkTriplePointRegion = triplePointRegion.grow(-2*(Point::Ones() - X - Y));

    for (auto iter : layout)
    {
        Box patchBox = layout[iter];
        Box localTriplePointRegion = shrunkTriplePointRegion & patchBox.adjacent(X+Y, ghostSize);
        for (auto center : localTriplePointRegion)
        {
            auto footprint = interpLayout.footprint(center, X+Y,iter);
            std::set<Point> footprintSln;
            Point corner = (X+Y)*(domainSize-1) + center*(Point::Ones()-X-Y);
            int dist = ((center-corner)*(X+Y)).abs().max();
            Point s0 = corner + X*dist - center;
            Point s1 = corner + Y*dist - center;

            for (auto fi : baseFootprint)
            {
                Point f0 = fi + s0;
                if (!triplePointRegion.containsPoint(center + f0))
                {
                    footprintSln.insert(f0);
                }
                Point f1 = fi + s1;
                if (!triplePointRegion.containsPoint(center + f1))
                {
                    footprintSln.insert(f1);
                }
            }
            for (auto fi : footprint)
            {
                EXPECT_TRUE(footprintSln.count(fi) == 1);
            }
        }
    }
}

#if DIM == 2
#if PR_VERBOSE > 0
TEST(MBInterpLayout, Visualization) {
    int boxSize = 8;
    int domainSize = boxSize;
    int ghostSize = 3;
    int numBlocks = 3;

    auto layout = testLayout(boxSize, domainSize, numBlocks);
    layout.print();
    MBInterpLayout interpLayout(layout);

    for (auto pi : Box::Kernel(2))
    {
        if (pi.sumAbs() <= 2) { interpLayout.addPoint(pi); }
    }
    MBLevelBoxData<int, 1, HOST> domainData(layout, Point::Zeros());
    domainData.setVal(-1);
    HDF5Handler h5;
    h5.writeLevel(domainData.getBlock(0), "DOMAIN");

    int plotIter = 0;
    for (auto iter : layout)
    {
        if (layout.block(iter) != 0)
        {
            continue;
        }
        Box patchBox = layout[iter];
        for (auto dir : Box::Kernel(1))
        {
            if (dir == Point::Zeros())
            {
                continue;
            }
            Box boundBox = patchBox.adjacent(dir * ghostSize);
            for (auto center : boundBox)
            {
                if (layout.domain().isPointInInterior(center, 0)) { continue; }
                if (layout.domain().isPointInDomainBoundary(center, 0)) { continue; }
                auto footprintData = printInterpLayout(interpLayout, center, dir, iter);

                h5.writePatch(footprintData, "FOOTPRINT_%i", plotIter);
                plotIter++;

            }
        }
    }
}
#endif
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
