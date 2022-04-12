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
    kern &= cube;
    EXPECT_EQ(kern==Box(Point::Ones(4)),true);
    Box grow = cube & Point::Basis(1,-2);
    EXPECT_EQ(grow.low()==Point(0,-2,0,0),true);
    cube&=grow.low();
    EXPECT_EQ(cube.low()==Point(0,-2,0,0),true);
    EXPECT_EQ(kern==Box(kern.low(),kern.high()),true);
    EXPECT_EQ(kern!=Box(kern.high(),kern.low()),true);
    EXPECT_EQ(cube<Box(Point::Ones()),true);
    EXPECT_EQ(kern%twos==twos,true);
    EXPECT_EQ(kern.mod(kern.low()-Point::Ones())==kern.high(),true);
    EXPECT_EQ(kern.mod(kern.high()+Point::Ones())==kern.low(),true);
    Box B1 = Box::Cube(2);
    Box B2 = B1.shift(0,1);
    Box B3 = B1.shift(1,-1);
    EXPECT_EQ(B2==Box(Point(1,0),Point(2,1)),true);
    EXPECT_EQ(B3==Box(Point(0,-1),Point(1,0)),true);
    B2 = B1.shift(Point({2,-3}));
    EXPECT_EQ(B2==Box(Point({2,-3}),Point({3,-2})),true);
    EXPECT_EQ(B2.toOrigin().low()==Point::Basis(0,0),true);
    B2 = B1.grow(3);
    EXPECT_EQ(B2==Box(Point::Ones(-3),Point::Ones(4)),true);
    EXPECT_EQ(B2.grow(-2)==Box(Point::Ones(-1),Point::Ones(2)),true);
    B1 = Box::Cube(4).grow(Point({-1,1}));
    EXPECT_EQ(B1==Box(Point({1,-1}),Point({2,4})),true);
    B1 = Box::Cube(4).grow(0,2);
    EXPECT_EQ(B1==Box(Point({-2,0}),Point(5,3)),true);
    B1 = Box::Cube(2);
    B2 = B1.grow(1,Side::Lo,2);
    EXPECT_EQ(B2==Box(Point::Basis(1,-2),Point::Ones()),true);
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
