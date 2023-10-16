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
