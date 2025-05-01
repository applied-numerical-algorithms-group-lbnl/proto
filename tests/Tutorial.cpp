#include <gtest/gtest.h>
#include "Proto.H"

using namespace Proto;

namespace {
    template<typename T, MemType MEM>
    PROTO_KERNEL_START
    void f_iotaTestF(Point& a_pt, Var<T,DIM,MEM>& a_X, T a_dx)
    {
        for (int dir = 0; dir < DIM; dir++)
        {
            a_X(dir) = a_pt[dir]*a_dx;
        }
    }
    PROTO_KERNEL_END(f_iotaTestF, f_iotaTest);
}

TEST(Tutorial, PointBox) {
    
    // This is a Point. Whe usually write them in pseudocode like this: (a,b,c,...);
    Point z = Point::Zeros();
    // Under the hood they are basically just a std::array<int, DIM>
    // DIM is a compile defined macro which defines the spatial dimensions of the build.
    
    // You can also build a Point like this:
    Point p1{1,2,3,4,5,6};
    // if DIM == 2, this Point is (1,2). if DIM == 3, this Point is (1,2,3)
    // this syntax prevents having a bunch if #if DIM {...} #endif code laying around
    // the above code will compile so long as DIM is less than the number of arguments

    // Here are some other useful ways to make Points:
    Point x = Point::X();       // (1,0,0,...)
    Point y = Point::Y();       // (0,1,0,...)
    Point b0 = Point::Basis(0); // Same as X
    Point u = Point::Ones();    // (1,1,1,...)

    EXPECT_EQ(x, Point(1,0,0,0,0,0));
    EXPECT_EQ(y, Point(0,1,0,0,0,0));
    EXPECT_EQ(x, b0);
    EXPECT_EQ(u, Point(1,1,1,1,1,1));

    // You can also index into them as if they were arrays (they are)
    EXPECT_EQ(x[0], 1);

    // Points support most arithmetic with other points and with scalars:
    Point uu = 2*u;
    Point xx = x + b0;

    EXPECT_EQ(uu, Point(2,2,2,2,2,2));
    EXPECT_EQ(xx, Point(2,0,0,0,0,0));

    // This is a Box. We usually write Boxes like [(i,j,k,...), (x,y,z,...)]
    // Boxes are DIM dimensional logically rectangular collections of Points
    // A Box is built by specifying its low corner and its high corner.
    // This Box is [(0,0,0,...), (1,1,1,...)]
    Box B(Point::Zeros(), Point::Ones());

    // IMPORTANT: The high corner of a Box is INCLUDED in its domain.
    // Explicitly, B above contains all Points which satisfy (0 <= p0 <= 1, 0 <= p1 <= 1, ...)
    EXPECT_TRUE(B.containsPoint(Point::Zeros()));
    EXPECT_TRUE(B.containsPoint(Point::Ones()));
    EXPECT_TRUE(B.containsPoint(Point::X()));
    EXPECT_TRUE(B.containsPoint(Point::Y()));

    // There are many many functions for mutating an existing Box:
    Box B1 = B.grow(1);                 // Add 1 to each range on the upper and lower sides
    Box B2 = B.extrude(Point::X(),1);   // Add 1 layer of Points in the positive X direction

    EXPECT_EQ(B1, Box(-u, 2*Point::Ones()));
    EXPECT_EQ(B2, Box(Point::Zeros(), Point::Ones() + Point::X()));

    // This is a BoxData. It contains N components of datatype double and it is allocated on the HOST (e.g. not a GPU)
    // It's domain indices are all of the points in Box B. 
    constexpr unsigned int N = DIM;
    BoxData<double, N, HOST> D(B);

    // We can initialize the data in a BoxData to a constant, or by indexing each Point one by one:
    D.setVal(0);
    for (auto pi : B)
    {
        for (int ii = 0; ii < N; ii++)
        {
            D(pi, ii) = pi.sum()*ii;
        }
    }

    // BoxData also usually support the kinds of arithmetic that you might expect (elementwise multiplication, scalar multiplication, etc)

    // We can also initialize a BoxData using forall. the f_iotaTest function is defined at the top of the file.
    double dx = 1.0/32;
    forallInPlace_p(f_iotaTest, D, dx);

    // We can also generate a new BoxData using forall
    auto D2 = forall_p<double, DIM>(f_iotaTest, B1, dx);
    // Note how the BoxData argument was replaced with a Box. That Box is the domain of the new BoxData:

    EXPECT_EQ(D2.box(), B1);
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
