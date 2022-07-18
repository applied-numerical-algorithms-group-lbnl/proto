#include <gtest/gtest.h>
#include "Proto.H"

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
    domain.close();
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
    hostData.initialize(f_MBPointID);
    for (auto iter : layout)
    {
        unsigned int block = layout.block(iter);
        Box box = layout[iter];
        BoxData<int, NCOMP, HOST>& hostData_i = hostData[iter];
        BoxData<int, NCOMP, HOST> slnData_i(box);
        forallInPlace_p(f_MBPointID, slnData_i, block);
        slnData_i -= hostData_i;
        EXPECT_TRUE(slnData_i.sum() == 0);
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
