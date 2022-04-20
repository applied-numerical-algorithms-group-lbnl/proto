#include <gtest/gtest.h>
#include <cmath>
#include "Proto.H"

using namespace Proto;
using namespace std;

TEST(BoxData, DefaultConstructor) {
    BoxData<double,2,MEMTYPE_DEFAULT,3> BD;
    EXPECT_TRUE(BD.box()==Box(Point::Zeros()));
    EXPECT_TRUE(BD.size()==0);
}

TEST(BoxData, BoxConstructor) {
    Box B = Box(Point(1,2,3,4,5,6,7));
    BoxData<int,3,MEMTYPE_DEFAULT,4,5> BD(B);
    EXPECT_TRUE(BD.box()==B);
    EXPECT_EQ(BD.size(),B.size()*3*4*5);
}

TEST(BoxData, Initializer) {
    Box B = Box(Point(1,2,3,4,5,6,7));
    BoxData<int,3,MEMTYPE_DEFAULT,4,5> BD(B,1337);
    EXPECT_EQ(BD.max(),1337);
    EXPECT_EQ(BD.min(),1337);
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
