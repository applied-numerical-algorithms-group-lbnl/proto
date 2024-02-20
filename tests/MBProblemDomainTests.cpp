#include <gtest/gtest.h>
#include "Proto.H"

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

TEST(MBProblemDomain, ConvertSimple)
{
    int boxSize = 4;
    for (auto dir : Box::Kernel(1))
    {
        if (dir.codim() != 1) { continue; }
        for (int axis = 0; axis < DIM; axis++)
        {
            if (dir[axis] != 0) { continue; }
            CoordPermutation R;
            for (int rot = 0; rot < 4; rot++)
            {
                MBProblemDomain domain(2);
                domain.defineBoundary(0,1,dir,R);
                domain.defineDomain(0,Point::Ones(boxSize));
                domain.defineDomain(1,Point::Ones(boxSize));
                const auto& graph = domain.graph();
                Point revDir = graph.reverseDir(0,1,dir);
                Point revArc = graph.reverseArc(0,1,dir);
                EXPECT_EQ(R(dir), -revDir);

                Box B0 = Box::Cube(boxSize).adjacent(dir);
                Box B1 = Box::Cube(boxSize).edge(revDir);
                Box B10 = domain.convert(B0,0,1);
                EXPECT_EQ(B10, B1);
                
                R = R*CoordPermutation::ccw(axis);
            }
        }
    }
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
// This test is turned off because Proto doesn't actually catch
// exceptions
#if 0
TEST(MBProblemDomain, Close) {
    CoordPermutation I;
    
    MBProblemDomain D0(2);
    D0.defineBoundary(0,1,0,Side::Hi,I);
    D0.defineDomain(0,Point::Ones(8));
    D0.defineDomain(1,Point::Ones(16));
    D0.close();

    MBProblemDomain D1(3);
    D1.defineBoundary(0,1,0,Side::Hi,I);
    D1.defineBoundary(1,2,1,Side::Hi,I);
    D1.defineDomain(0, Point{ 8, 8, 8, 8, 8, 8});
    D1.defineDomain(1, Point{16, 8, 8, 8, 8, 8});
    D1.defineDomain(2, Point{16,16,16,16,16,16});
    D1.close();
}
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
