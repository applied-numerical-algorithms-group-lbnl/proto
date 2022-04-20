#include <gtest/gtest.h>
#include <cmath>
#include "Proto.H"

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
    Box B0 = Box::Cube(4).shift(Point::Ones());
    B1 = B0.extrude(Point::Ones(), 2);
    B2 = B0.extrude(Point::Basis(0,-1),3);
    B3 = B0.extrude(Point(-1, 1, 0));
    EXPECT_EQ(B1==Box(Point::Ones(),Point::Ones(6)),true);
    EXPECT_EQ(B2==Box(Point(-2,1,1),Point::Ones(4)),true);
    EXPECT_EQ(B3==Box(Point(0,1,1),Point(4,5,4)),true);
    B1 = B0.extrude(1,3,false);
    B2 = B0.extrude(1,3,true);
    EXPECT_EQ(B1==Box(Point(1,-2,1),Point(4,4,4)),true);
    EXPECT_EQ(B2==Box(Point(1,1,1),Point(4,7,4)),true);
    EXPECT_EQ(B1.grow(1,Side::Hi,2)==B1.growHi(1,2),true);
    EXPECT_EQ(B2.grow(0,Side::Lo,1)==B2.growLo(0,1),true);
    B0 = Box::Cube(4);
    B1 = B0.coarsen(2);
    EXPECT_EQ(B1==Box::Cube(2),true);
    B2 = Box::Cube(3).shift(Point::Ones(2));
    B3 = B2.coarsen(2);
    EXPECT_EQ(B3==Box(Point::Ones(2)).shift(Point::Ones()),true);
    B1 = B0.coarsen(Point({1,2}));
    EXPECT_EQ(B1==Box(Point::Zeros(),Point(3,1,3)),true);
    EXPECT_EQ(B0.coarsenable(2),true);
    EXPECT_EQ(B0.coarsenable(Point(3,1,3)),false);
    EXPECT_EQ(B0.taperCoarsen(Point::Ones(2))==B0.coarsen(2),true);
    B1 = B0.shift(Point::Ones()).taperCoarsen(Point::Ones(2));
    EXPECT_EQ(B1==Box(B1.low(),Point::Ones(2)),true);
    B0 = Box::Cube(2);
    EXPECT_EQ(B0.refine(Point({1,2}))==Box(Point::Zeros(),Point(1,3)),true);
    B2 = Box::Cube(2).shift(Point::Ones());
    EXPECT_EQ(B2.refine(2)==Box(Point::Ones(2),Point::Ones(5)),true);
    B0 = Box::Cube(4).shift(Point::Ones()); //[(1,1), (4,4)]
    B1 = B0.edge(Point::Basis(0));          //[(4,1), (4,4)]
    B2 = B0.edge(Point::Ones(), 2);         //[(3,3), (4,4)]
    B3 = B0.edge(Point::Basis(1,-1), 2);    //[(1,1), (4,2)]
    EXPECT_EQ(B1==Box(Point(4,1,1),Point::Ones(4)),true);
    EXPECT_EQ(B2==Box(Point::Ones(3),Point::Ones(4)),true);
    EXPECT_EQ(B3==Box(Point::Ones(),Point(4,2,4)),true);
    B0 = Box::Cube(4).shift(Point::Ones()); //[(1,1), (4,4)]
    B1 = B0.face(0, Side::Hi);              //[(4,1), (4,4)]
    B2 = B0.face(1, Side::Lo, 2);           //[(1,1), (4,2)]
    B3 = B0.face(0, Side::Lo, 5);           //[(1,1), (4,4)]
    EXPECT_EQ(B1==Box(Point(4,1,1),Point::Ones(4)),true);
    EXPECT_EQ(B2==Box(Point::Ones(),Point(4,2,4)),true);
    EXPECT_EQ(B0==B3,true);
    EXPECT_EQ(B1==B0.flatten(0,true),true);
    EXPECT_EQ(B0.face(1,Side::Lo)==B0.flatten(1),true);
    B0 = Box::Cube(8);  // [(0,0) ,  (7,7)]
    B1 = B0.adjacent(Point(1,0,0) , 2); // [(8,0) ,  (9,7)]
    B2 = B0.adjacent(Point(0,-1,0), 2); // [(0,-2), (7,-1)]
    B3 = B0.adjacent(Point(-1,1,0), 2); // [(-2,8), (-1,9)]
    Box B4 = B0.adjacent(Point(1,0,0));     // [(8,0) , (15,7)]
    EXPECT_EQ(B1==Box(Point::Basis(0,8),Point(9,7,7)),true);
    EXPECT_EQ(B2==Box(Point::Basis(1,-2),Point(7,-1,7)),true);
    EXPECT_EQ(B3==Box(Point(-2,8,0),Point(-1,9,7)),true);
    EXPECT_EQ(B4==Box(Point::Basis(0,8),Point(15,7,7)),true);
    EXPECT_EQ(B1==B0.adjacent(0, Side::Hi, 2),true); // [(8, 0), (9, 7)]
    EXPECT_EQ(B2==B0.adjacent(1, Side::Lo, 2),true); // [(0,-2), (7,-1)]
    EXPECT_EQ(B4==B0.adjacent(0, Side::Hi),true);    // [(8, 0), (15,7)]
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
