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
  
    /*
    int numTests = 0;
    for (int bi = 0; bi < XPOINT_SIZE; bi++)
    {
        auto blockLayout = layout.layout(bi);
        Box blockDomainBox = blockLayout.domain().box();
        Box K = Box::Kernel(1);
        for (auto dir : K)
        {
            if (dir == Point::Zeros()) { continue; }
            auto bounds = hostData.bounds(bi, dir);
            Box localBoundBox = blockDomainBox.adjacent(dir,ghostSize);
            if (bounds.size() > 0)
            {
                EXPECT_TRUE(dir == nx || dir == ny || dir == nx+ny);
                for (auto b : bounds)
                {
                    EXPECT_EQ(b.localBlock, bi);
                    EXPECT_NE(b.adjBlock, bi);
                    //std::cout << b.localData->box() << " == " << localBoundBox << std::endl;
                    EXPECT_EQ(b.localData->box(), localBoundBox);
                }
                unsigned int xBlock = (bi+1) % XPOINT_SIZE;
                unsigned int yBlock = (bi-1+XPOINT_SIZE) % XPOINT_SIZE;
                if (dir == nx)
                {
                    Box adjBoundBox = domain.convert(localBoundBox, bi, xBlock);
                    EXPECT_EQ(bounds.size(), 1);
                    EXPECT_EQ(bounds[0].adjBlock, xBlock);
                    //std::cout << bounds[0].adjData->box() << " == " << adjBoundBox << std::endl;
                    EXPECT_EQ(bounds[0].adjData->box(), adjBoundBox);
                } else if (dir == ny)
                {
                    Box adjBoundBox = domain.convert(localBoundBox, bi, yBlock);
                    EXPECT_EQ(bounds.size(), 1);
                    EXPECT_EQ(bounds[0].adjBlock, yBlock);
                    //std::cout << bounds[0].adjData->box() << " == " << adjBoundBox << std::endl;
                    EXPECT_EQ(bounds[0].adjData->box(), adjBoundBox);
                } else if (dir == nx+ny)
                {
                    EXPECT_EQ(bounds.size(), max(XPOINT_SIZE-3, 0));
                    for (auto b : bounds)
                    {
                        EXPECT_NE(b.adjBlock, xBlock);
                        EXPECT_NE(b.adjBlock, yBlock);
                        Box cornerBox = Box::Cube(domainSize).edge(nx+ny, ghostSize);
                        //std::cout << b.adjData->box() << " == " << cornerBox << std::endl;
                        EXPECT_EQ(b.adjData->box(), cornerBox);
                    }
                } 
            }
        }
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
