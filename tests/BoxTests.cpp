#include <gtest/gtest.h>
#include <cmath>
#include "Proto.H"

#define _USE_MATH_DEFINES

using namespace Proto;
using namespace std;

TEST(Base, Box) {
    Box cube = Box::Cube(4);
    EXPECT_EQ(cube.empty(),false);
    EXPECT_EQ(cube.size(),int(pow(4,DIM)));
    EXPECT_EQ(cube.contains(Box(Point::Ones(4))),true);
    for (int i=0; i<DIM; i++)
        EXPECT_EQ(cube.size(i),4);
    EXPECT_EQ(cube.onBoundary(Point(0,2,3)),true);
    Box kern = Box::Kernel(3);
    EXPECT_EQ(kern.low()==Point(-3,-3,-3,-3),true);
    EXPECT_EQ(kern.high()==Point(3,3,3,3),true);
    Point p(-3,-2,-1,0);
    size_t idx=0;
    for (size_t i=0; i<DIM; i++)
        idx += size_t(pow(kern.size(i),i)) * (p[i]-kern.low()[i]);
    EXPECT_EQ(kern.index(p),idx);
    EXPECT_EQ(kern[idx]==p,true);
    Point low = kern.low();
    unsigned int z = 0;
#if DIM == 3 
    z = p[2]-low[2];
#endif
    EXPECT_EQ(kern(p[0]-low[0],p[1]-low[1],z)==p,true);
    Point twos = Point::Ones(2);
    EXPECT_EQ((cube&cube.shift(twos)) == Box(twos,Point::Ones(3)),true);
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
