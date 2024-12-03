#include <gtest/gtest.h>
#include "Proto.H"

using namespace Proto;

TEST(Point, VectorSpace) {
    Point p0 = Point::Basis(0);
    Point p1 = Point::Basis(1,-1);
    Point p2 = p0 + p1;

    auto par0 = p0.parallelUnit();
    auto par1 = p1.parallelUnit();
    auto par2 = p2.parallelUnit();
    
    auto prp0 = p0.perpUnit();
    auto prp1 = p1.perpUnit();
    auto prp2 = p2.perpUnit();

    for (int dir = 0; dir < 2; dir++)
    {
        if (dir == 0)
        {
            EXPECT_EQ(par0[dir], Point::Basis(0));
            EXPECT_EQ(par1[dir], Point::Zeros());
            EXPECT_EQ(par2[dir], Point::Basis(0));
            
            EXPECT_EQ(prp0[dir], Point::Zeros());
            EXPECT_EQ(prp1[dir], Point::Basis(0));
            EXPECT_EQ(prp2[dir], Point::Zeros());
        } else if (dir == 1)
        {
            EXPECT_EQ(par0[dir], Point::Zeros());
            EXPECT_EQ(par1[dir], Point::Basis(1));
            EXPECT_EQ(par2[dir], Point::Basis(1));
            
            EXPECT_EQ(prp0[dir], Point::Basis(1));
            EXPECT_EQ(prp1[dir], Point::Zeros());
            EXPECT_EQ(prp2[dir], Point::Zeros());
        } else {
            EXPECT_EQ(par0[dir], Point::Basis(dir));
            EXPECT_EQ(par1[dir], Point::Basis(dir));
            EXPECT_EQ(par2[dir], Point::Basis(dir));
        }
    }
}

TEST(Point, BinaryReduction) {
    Point p1{ 1, 2, 3, 4,  5,  6};
    Point p2{-2, 4,-6, 8,-10, 12};
    Point pmax_0{1,4,3,8,5,12};
    Point pmin_0{-2,2,-6,4,-10,6};
    Point pabsMax_0{2,4,6,8,10,12};

    auto pmax = maxPoint(p1,p2);
    auto pmin = minPoint(p1,p2);
    auto pabsMax = absMaxPoint(p1,p2);

    EXPECT_EQ(pmax, pmax_0);
    EXPECT_EQ(pmin, pmin_0);
    EXPECT_EQ(pabsMax, pabsMax_0);

}

TEST(Point, DirectionsOfCodim)
{
    #if DIM == 2
        std::set<Point> codim1Soln;
        std::set<Point> codim2Soln;

        codim1Soln.insert(Point(+1,0));
        codim1Soln.insert(Point(-1,0));
        codim1Soln.insert(Point(0,+1));
        codim1Soln.insert(Point(0,-1));

        codim2Soln.insert(Point(+1,+1));
        codim2Soln.insert(Point(+1,-1));
        codim2Soln.insert(Point(-1,+1));
        codim2Soln.insert(Point(-1,-1));

        auto codim1Test = Point::DirectionsOfCodim(1);
        EXPECT_EQ(codim1Soln.size(), codim1Test.size());
        for (auto pi : codim1Soln)
        {
            EXPECT_TRUE(codim1Test.find(pi) != codim1Test.end());
        }

        auto codim2Test = Point::DirectionsOfCodim(2);
        EXPECT_EQ(codim2Soln.size(), codim2Test.size());
        for (auto pi : codim2Soln)
        {
            EXPECT_TRUE(codim2Test.find(pi) != codim2Test.end());
        }
    #endif
    #if DIM == 3
        std::set<Point> codim1Soln;
        std::set<Point> codim2Soln;
        std::set<Point> codim3Soln;

        codim1Soln.insert(Point(+1,0,0));
        codim1Soln.insert(Point(-1,0,0));
        codim1Soln.insert(Point(0,+1,0));
        codim1Soln.insert(Point(0,-1,0));
        codim1Soln.insert(Point(0,0,+1));
        codim1Soln.insert(Point(0,0,-1));

        codim2Soln.insert(Point(+1,+1,0));
        codim2Soln.insert(Point(+1,-1,0));
        codim2Soln.insert(Point(-1,+1,0));
        codim2Soln.insert(Point(-1,-1,0));
        codim2Soln.insert(Point(+1,0,+1));
        codim2Soln.insert(Point(+1,0,-1));
        codim2Soln.insert(Point(-1,0,+1));
        codim2Soln.insert(Point(-1,0,-1));
        codim2Soln.insert(Point(0,+1,+1));
        codim2Soln.insert(Point(0,-1,-1));
        codim2Soln.insert(Point(0,-1,+1));
        codim2Soln.insert(Point(0,+1,-1));

        codim3Soln.insert(Point(+1,+1,+1));
        codim3Soln.insert(Point(+1,-1,+1));
        codim3Soln.insert(Point(-1,+1,+1));
        codim3Soln.insert(Point(-1,-1,+1));
        codim3Soln.insert(Point(+1,+1,-1));
        codim3Soln.insert(Point(+1,-1,-1));
        codim3Soln.insert(Point(-1,+1,-1));
        codim3Soln.insert(Point(-1,-1,-1));

        auto codim1Test = Point::DirectionsOfCodim(1);
        EXPECT_EQ(codim1Soln.size(), codim1Test.size());
        for (auto pi : codim1Soln)
        {
            EXPECT_TRUE(codim1Test.find(pi) != codim1Test.end());
        }

        auto codim2Test = Point::DirectionsOfCodim(2);
        EXPECT_EQ(codim2Soln.size(), codim2Test.size());
        for (auto pi : codim2Soln)
        {
            EXPECT_TRUE(codim2Test.find(pi) != codim2Test.end());
        }

        auto codim3Test = Point::DirectionsOfCodim(3);
        EXPECT_EQ(codim3Soln.size(), codim3Test.size());
        for (auto pi : codim3Soln)
        {
            EXPECT_TRUE(codim3Test.find(pi) != codim3Test.end());
        }
    #endif
}

TEST(Point, Directions)
{
    auto dirs = Point::Directions();
    std::set<Point> correctDirs;
    for (auto p : Box::Kernel(1))
    {
        if (p == Point::Zeros()) { continue; }
        correctDirs.insert(p);
    }
    EXPECT_EQ(correctDirs.size(), dirs.size());
    for (auto dir : dirs)
    {
        EXPECT_NE(correctDirs.find(dir), correctDirs.end());
    }
}

TEST(Point, NonZeroIndex)
{
    Point p0 = Point::Zeros();
    Point p1 = Point::Basis(0,-1);
    Point p2 = Point::Basis(DIM-1);
    Point p3 = Point::Ones();

    EXPECT_EQ(p0.firstNonZeroIndex(), DIM);
    EXPECT_EQ(p1.firstNonZeroIndex(), 0);
    EXPECT_EQ(p2.firstNonZeroIndex(), DIM-1);
    EXPECT_EQ(p3.firstNonZeroIndex(), 0);

    EXPECT_EQ(p0.lastNonZeroIndex(), -1);
    EXPECT_EQ(p1.lastNonZeroIndex(), 0);
    EXPECT_EQ(p2.lastNonZeroIndex(), DIM-1);
    EXPECT_EQ(p3.lastNonZeroIndex(), DIM-1);
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
