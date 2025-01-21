#include <gtest/gtest.h>
#include "Proto.H"
#include "TestFunctions.H"
using namespace Proto;

namespace {
    MBDisjointBoxLayout testLayout(int boxSize, int domainSize, int numBlocks)
    {
        auto domain = buildXPoint(domainSize, numBlocks);
        std::vector<MBPoint> patches;
        std::vector<Point> boxSizes;
        for (BlockIndex bi = 0; bi < numBlocks; bi++)
        {
            for (auto pi : Box::Cube(domainSize / boxSize))
            {
                if (pi != Point::X())
                {
                    patches.push_back(MBPoint(pi, bi));
                }
            }
            boxSizes.push_back(Point::Ones(boxSize));
        }
        return MBDisjointBoxLayout(domain, patches, boxSizes);
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

TEST(MBInterpLayout, BlockBoundary) {
    int boxSize = 8;
    int domainSize = boxSize*2;
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
#if PR_VERBOSE > 0
    HDF5Handler h5;
    h5.writeLevel(domainData.getBlock(0), "DOMAIN");
#endif

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
#if PR_VERBOSE > 0
                h5.writePatch(footprintData, "FOOTPRINT_%i", plotIter);
                plotIter++;
#endif
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
