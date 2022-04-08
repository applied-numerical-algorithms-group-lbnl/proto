#include <gtest/gtest.h>
#include <cmath>
#include "Proto.H"

#define _USE_MATH_DEFINES

using namespace Proto;
using namespace std;

TEST(Box, Cube) {
    Box cube = Box::Cube(4);
    EXPECT_EQ(cube.empty(),false);
    EXPECT_EQ(cube.size(),int(pow(4,DIM)));
    EXPECT_EQ(cube.contains(Box(Point::Ones(4))),true);
    for (int i=0; i<DIM; i++)
        EXPECT_EQ(cube.size(i),4);
    EXPECT_EQ(cube.onBoundary(Point(0,2,3)),true);
}

TEST(Box, Kernel) {
    Box kern = Box::Kernel(3);
    EXPECT_EQ(kern.low()==Point(-3,-3,-3,-3),true);
    EXPECT_EQ(kern.high()==Point(3,3,3,3),true);
    Point p(-3,-2,-1,0);
    int idx=0;
    for (int i=0; i<DIM; i++)
        idx += int(pow(kern.size(i),i)) * (p[i]-kern.low()[i]);
    EXPECT_EQ(kern.index(p),idx);
    EXPECT_EQ(kern[idx]==p,true);
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
