#include <gtest/gtest.h>
#include "ProtoMMB.H"

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

TEST(MBProblemDomain, Convert) {
    int domainSize = 64;
    auto domain = buildXPoint(domainSize);
    Point x = Point::Basis(0);
    Point y = Point::Basis(1);
    Point origin = Point::Zeros();
    Box domainBox = Box::Cube(domainSize);
    Box xAdj = domainBox.adjacent(x, 1);
    Box yAdj = domainBox.adjacent(y, 1);
    Box xEdge = domainBox.edge(x, 1);
    Box yEdge = domainBox.edge(y, 1);
    std::cout << "xEdge: " << xEdge << " | yEdge: " << yEdge << std::endl;
    for (unsigned int bi = 0; bi < domain.size(); bi++)
    {
        unsigned int bx = (bi + 1) % XPOINT_SIZE;
        unsigned int by = (bi + XPOINT_SIZE - 1) % XPOINT_SIZE;
        Point sx = x*domainSize*2;
        Point sy = y*domainSize*2;

        EXPECT_EQ(domain.convert(origin, bi, bx), sy);
        EXPECT_EQ(domain.convert(origin, bi, by), sx);
    
        EXPECT_EQ(domain.convert(xAdj, bi, bx), yEdge);
        EXPECT_EQ(domain.convert(yAdj, bi, by), xEdge);
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
