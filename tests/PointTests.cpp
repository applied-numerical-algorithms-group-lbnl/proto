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

    for (int dir = 0; dir < DIM; dir++)
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
