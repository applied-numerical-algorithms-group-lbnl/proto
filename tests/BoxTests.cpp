#include <gtest/gtest.h>
#include <cmath>
#include "Proto.H"

using namespace Proto;
using namespace std;

// TODO: Split this up...
TEST(Box, Base) {
    if (DIM > 6)
    {
        MayDay<void>::Warning("Some tests may not pass with DIM > 6");
    }
    Box cube = Box::Cube(4);
    EXPECT_FALSE(cube.empty());
    EXPECT_EQ(cube.size(),int(pow(4,DIM)));
    EXPECT_TRUE(cube.containsBox(Box(Point::Ones(4))));
    for (int i=0; i<DIM; i++)
        EXPECT_EQ(cube.size(i),4);
    EXPECT_TRUE(cube.onBoundary(Point(0,2,3)));
    Box kern = Box::Kernel(3);
    EXPECT_EQ(kern.low(),Point(-3,-3,-3,-3));
    EXPECT_EQ(kern.high(),Point(3,3,3,3));
    Point p(-3,-2,-1,0);
    size_t idx=0;
    for (size_t i=0; i<DIM; i++)
        idx += size_t(pow(kern.size(i),i)) * (p[i]-kern.low()[i]);
    EXPECT_EQ(kern.index(p),idx);
    EXPECT_EQ(kern[idx],p);
    Point low = kern.low();
    unsigned int z = 0;
#if DIM == 3 
    z = p[2]-low[2];
#endif
    EXPECT_EQ(kern(p[0]-low[0],p[1]-low[1],z),p);
    Point twos = Point::Ones(2);
    EXPECT_EQ((cube&cube.shift(twos)) , Box(twos,Point::Ones(3)));
    kern &= cube;
    EXPECT_EQ(kern,Box(Point::Ones(4)));
    Box grow = cube & Point::Basis(1,-2);
    EXPECT_EQ(grow.low(),Point(0,-2,0,0));
    cube&=grow.low();
    EXPECT_EQ(cube.low(),Point(0,-2,0,0));
    EXPECT_EQ(kern,Box(kern.low(),kern.high()));
    EXPECT_NE(kern,Box(kern.high(),kern.low()));
    EXPECT_LT(cube,Box(Point::Ones()));
    EXPECT_EQ(kern%twos,twos);
    EXPECT_EQ(kern.mod(kern.low()-Point::Ones()),kern.high());
    EXPECT_EQ(kern.mod(kern.high()+Point::Ones()),kern.low());
    Box B1 = Box::Cube(2);
    Box B2 = B1.shift(0,1);
    Box B3 = B1.shift(1,-1);
    EXPECT_EQ(B2,Box(Point(1,0,0,0,0,0),Point(2,1,1,1,1,1)));
    EXPECT_EQ(B3,Box(Point(0,-1,0,0,0,0),Point(1,0,1,1,1,1)));
    B2 = B1.shift(Point(2,-3,0,0,0,0));
    EXPECT_EQ(B2,Box(Point(2,-3,0,0,0,0),Point(3,-2,1,1,1,1)));
    EXPECT_EQ(B2.toOrigin().low(),Point::Basis(0,0));
    B2 = B1.grow(3);
    EXPECT_EQ(B2,Box(Point::Ones(-3),Point::Ones(4)));
    EXPECT_EQ(B2.grow(-2),Box(Point::Ones(-1),Point::Ones(2)));
    B1 = Box::Cube(4).grow(Point(-1,1,2,3,4,5));
    EXPECT_EQ(B1,Box(Point(1,-1,-2,-3,-4,-5),Point(2,4,5,6,7,8)));
    B1 = Box::Cube(4).grow(0,2);
    EXPECT_EQ(B1,Box(Point(-2,0,0,0,0,0),Point(5,3,3,3,3,3)));
    B1 = Box::Cube(2);
    B2 = B1.grow(1,Side::Lo,2);
    EXPECT_EQ(B2,Box(Point::Basis(1,-2),Point::Ones()));
    Box B0 = Box::Cube(4).shift(Point::Ones());
    B1 = B0.extrude(Point::Ones(), 2);
    B2 = B0.extrude(Point::Basis(0,-1),3);
    B3 = B0.extrude(Point(-1, 1, 0));
    EXPECT_EQ(B1,Box(Point::Ones(),Point::Ones(6)));
    EXPECT_EQ(B2,Box(Point(-2,1,1),Point::Ones(4)));
    EXPECT_EQ(B3,Box(Point(0,1,1),Point(4,5,4)));
    B1 = B0.extrude(1,3,false);
    B2 = B0.extrude(1,3);
    EXPECT_EQ(B1,Box(Point(1,-2,1),Point(4,4,4)));
    EXPECT_EQ(B2,Box(Point(1,1,1),Point(4,7,4)));
    EXPECT_EQ(B1.grow(1,Side::Hi,2),B1.growHi(1,2));
    EXPECT_EQ(B2.grow(0,Side::Lo,1),B2.growLo(0,1));
    B0 = Box::Cube(4);
    B1 = B0.coarsen(2);
    EXPECT_EQ(B1,Box::Cube(2));
    B2 = Box::Cube(3).shift(Point::Ones(2));
    B3 = B2.coarsen(2);
    EXPECT_EQ(B3,Box(Point::Ones(2)).shift(Point::Ones()));
    B0 = Box::Cube(64);
    B1 = B0.coarsen(Point(1,2,4,8,16,32));
    EXPECT_EQ(B1,Box(Point::Zeros(),Point(63,31,15,7,3,1)));
    EXPECT_TRUE(B0.coarsenable(2));
    EXPECT_TRUE(B0.coarsenable(Point(1,2,4,8,16,32)));
    EXPECT_EQ(B0.taperCoarsen(Point::Ones(2)),B0.coarsen(2));
    B1 = B0.shift(Point::Ones()).taperCoarsen(Point::Ones(2));
    EXPECT_EQ(B1,Box(B1.low(),Point::Ones(32)));
    B0 = Box::Cube(2).shift(Point::Ones());
    B1 = B0.refine(Point(1,2,3,4,5,6));
    EXPECT_EQ(B1, Box(Point(1,2,3,4,5,6),Point(2,5,8,11,14,17)));
    B2 = Box::Cube(2).shift(Point::Ones());
    EXPECT_EQ(B2.refine(2),Box(Point::Ones(2),Point::Ones(5)));
    B0 = Box::Cube(4).shift(Point::Ones()); //[(1,1), (4,4)]
    B1 = B0.face(0, Side::Hi);              //[(4,1), (4,4)]
    B2 = B0.face(1, Side::Lo, 2);           //[(1,1), (4,2)]
    B3 = B0.face(0, Side::Lo, 5);           //[(1,1), (4,4)]
    EXPECT_EQ(B1,Box(Point(4,1,1),Point::Ones(4)));
    EXPECT_EQ(B2,Box(Point::Ones(),Point(4,2,4)));
    EXPECT_EQ(B0,B3);
    EXPECT_EQ(B1,B0.flatten(0,true));
    EXPECT_EQ(B0.face(1,Side::Lo),B0.flatten(1));
}

TEST(Box, Adjacent) {
    Point x = Point::Basis(0);
    Point y = Point::Basis(1);
    int D = 8;
    Box B0 = Box::Cube(D);        // [(0,0) ,  (7,7)]
    Point L = B0.low();
    Point H = B0.high();

    Box B1 = B0.adjacent(x , 2);  // [(8,0) ,  (9,7)]
    EXPECT_EQ(B1, Box(L + D*x,          H + 2*x));
    Box B2 = B0.adjacent(-y, 2);  // [(0,-2), (7,-1)]
    EXPECT_EQ(B2, Box(L - 2*y,          H - D*y));
    Box B3 = B0.adjacent(y-x, 2); // [(-2,8), (-1,9)]
    EXPECT_EQ(B3, Box(L - 2*x + D*y,  H - D*x + 2*y));
    Box B4 = B0.adjacent(x,-1);   // [(7,0) , (7,7)]
    //EXPECT_EQ(B4, Box(L + D*x,  H + D*x));
    EXPECT_EQ(B4, Box(L+(D-1)*x, H));
    Box B5 = B0.adjacent(x-(2*y));    // [(8,-2) , (15,-1)]
    EXPECT_EQ(B5, Box(L + D*x - 2*y,  H + x - D*y));
}

TEST(Box, WhichBoundaryContains)
{
    int ghostSize = 2;
    Box B0 = Box::Kernel(2);
    Box B1 = B0.grow(ghostSize);

    for (auto dir : Box::Kernel(1))
    {
        auto boundBox = B0.adjacent(dir, ghostSize);
        for (auto pi : boundBox)
        {
            auto boundDir = B0.whichBoundaryContains(pi);
            EXPECT_EQ(boundDir, dir);
        }
    }
}

TEST(Box, Edge) {
#if DIM > 3
    MayDay<void>::Warning("Box.Edge test was not designed to be run for DIM>3");
#endif

    Point x = Point::Basis(0);
    Point y = Point::Basis(1);
    int D = 8;
    Box B0 = Box::Cube(D).shift(Point::Ones());        // [(1,1) ,  (8,8)]
    Point L = B0.low();
    Point H = B0.high();

    Box B1 = B0.edge(Point::Basis(0));          //[(8,1,1), (8,8,8)]
    Box B2 = B0.edge(Point::Ones(), 2);         //[(7,7,7), (8,8,8)]
    Box B3 = B0.edge(x-y, 2);                   //[(7,1,1), (8,2,8)]
    EXPECT_EQ(B1,Box(L + x*(D-1), H));
    EXPECT_EQ(B2,Box(L + Point::Ones()*(D-2), H));
    EXPECT_EQ(B3,Box(L + x*(D-2), H - y*(D-2)));
}

TEST(Box, UnionOperator) {
    Box B0 = Box::Cube(8);
    Box B1 = Box::Cube(8).shift(Point::Ones());
    Box U0 = B0 + B1;
    EXPECT_EQ(U0, Box::Cube(9));

    Box B2 = Box::Cube(8).shift(-Point::Ones());
    B2 += B1;
    EXPECT_EQ(B2, Box::Cube(10).shift(-Point::Ones()));
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
