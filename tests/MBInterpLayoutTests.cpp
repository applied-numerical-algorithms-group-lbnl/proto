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
            // std::cout << "checking center " << center << std::endl;
            auto footprint = interpLayout.footprint(center, X, iter);
            for (auto fi : xBoundFootprint0)
            {
                // std::cout << "\tlooking for shift: " << fi << std::endl;
                EXPECT_TRUE(footprint.count(fi) == 1);
            }
        }
        for (auto center : xBound1 & localXBound)
        {
            // std::cout << "checking center " << center << std::endl;
            auto footprint = interpLayout.footprint(center, X, iter);
            for (auto fi : xBoundFootprint1)
            {
                // std::cout << "\tlooking for shift: " << fi << std::endl;
                EXPECT_TRUE(footprint.count(fi) == 1);
            }
        }
        auto localYBound = layout[iter].adjacent(Y, ghostSize);
        for (auto center : yBound0 & localYBound)
        {
            // std::cout << "checking center " << center << std::endl;
            auto footprint = interpLayout.footprint(center, Y, iter);
            for (auto fi : yBoundFootprint0)
            {
                // std::cout << "\tlooking for shift: " << fi << std::endl;
                EXPECT_TRUE(footprint.count(fi) == 1);
            }
        }
        for (auto center : yBound1 & localYBound)
        {
            // std::cout << "checking center " << center << std::endl;
            auto footprint = interpLayout.footprint(center, Y, iter);
            for (auto fi : yBoundFootprint1)
            {
                // std::cout << "\tlooking for shift: " << fi << std::endl;
                EXPECT_TRUE(footprint.count(fi) == 1);
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
