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

TEST(SUITE_NAME, TEST_NAME) {
    int domainSize = 64;
    int boxSize = 16;
    int ghostSize = 1;
    auto domain = buildXPoint(domainSize);
    Point boxSizeVect = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);

    MBLevelBoxData<int, NCOMP, HOST> hostData(layout, ghostSize);
   
    /*
    for (int bi = 0; bi < XPOINT_SIZE; bi++)
    {
        auto blockLayout = layout.layout(bi);
        Box blockDomainBox = blockLayout.domain().box();

        for (auto iter : layout)
        {
            auto localIndex = layout.blockIndex(iter);
            auto& patch = hostData[iter];
            Point patchID = blockLayout.point(localIndex);
            Box K = Box::Kernel(1);
            for (auto dir : K)
            {
                if (dir == Point::Zeros()) { continue; }
                Box edgeBox = blockDomainBox.edge(dir, 1);
                if (edgeBox.contains(patchID))
                {
                    Box localBox = patch.box().adjacent(dir,1);
                    auto bounds = domain.graph().boundaries(bi, dir);
                    for (auto bound : bounds)
                    {
                        auto adjBlock = bound.dstBlock;
                        auto adjBox = domain.convert(localBox, bi, adjBlock);

                    }
                }
            }

        }
    }
    hostData.initialize(f_MBPointID);
    for (auto iter : layout)
    {
        unsigned int block = layout.block(iter);
        Box box = layout[iter];
        BoxData<int, NCOMP, HOST>& hostData_i  = hostData[iter];
        BoxData<int, NCOMP, HOST> slnData_i(box);
        forallInPlace_p(f_MBPointID, slnData_i, block);
        slnData_i -= hostData_i;
        EXPECT_TRUE(slnData_i.sum() == 0);
    }
    */
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
