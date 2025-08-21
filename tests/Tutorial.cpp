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
    
    //=============================================================================================
    // POINT

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

    //=============================================================================================
    // BOX

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

    // We use the Cube function a lot. It creates a Box with low point at zero with specified side length
    auto B1 = Box::Cube(32);
    EXPECT_EQ(B1.low(), Point::Zeros());
    EXPECT_EQ(B1.high(), 31*Point::Ones()); //Note the Box includes (31, 31, 31, ...) but not 32

    // There are many many functions for mutating an existing Box:
    Box B2 = B.grow(1);                 // Add 1 to each range on the upper and lower sides
    Box B3 = B.extrude(Point::X(),1);   // Add 1 layer of Points in the positive X direction

    EXPECT_EQ(B2, Box(-u, 2*Point::Ones()));
    EXPECT_EQ(B3, Box(Point::Zeros(), Point::Ones() + Point::X()));

    // You can also intersect boxes:
    Box A1 = Box::Cube(32);
    Box A2 = Box::Cube(64);
    auto A3 = A1 & A2; 
    EXPECT_EQ(A1, A3);

    //=============================================================================================
    // BOXDATA

    // This is a BoxData. It contains N components of datatype double and it is allocated on the HOST (e.g. not a GPU)
    // It's domain indices are all of the points in Box B. 
    constexpr unsigned int N = DIM;
    BoxData<double, N, HOST> D(B);

    // We can initialize the data in a BoxData to a constant, or by indexing each Point one by one:
    D.setVal(0);
    for (auto pi : D.box())
    {
        for (int ii = 0; ii < N; ii++)
        {
            D(pi, ii) = pi.sum()*ii;
        }
    }

    // BoxData also usually support the kinds of arithmetic that you might expect (elementwise multiplication, scalar multiplication, etc)

    //=============================================================================================
    // FORALL

    // We can also initialize a BoxData using forall. the f_iotaTest function is defined at the top of the file.
    double dx = 1.0/32;
    forallInPlace_p(f_iotaTest, D, dx);

    // We can also generate a new BoxData using forall
    auto D2 = forall_p<double, DIM>(f_iotaTest, B1, dx);
    
    // Note how the BoxData argument was replaced with a Box. That Box is the domain of the new BoxData:
    EXPECT_EQ(D2.box(), B1);

    // You can also use forall with lambdas as follows.
    BoxData<double, 1> D3(Box::Cube(32), 1.0);
    BoxData<double, 1> D4(Box::Cube(32).shift(Point::Ones())); //shift adds the specified point to both low and high, translating the Box

    // Note the PROTO_LAMBDA macro
    auto D5 = forall_p<double>
    ([] PROTO_LAMBDA(Point& p, Var<double, 1>& v5, Var<double, 1>& v3, Var<double, 1>& v4)
        {
            v5(0) = v3(0) + v4(0); //summing the 0th (only) component
        },
    D3,D4);
    // NB: the "_p" suffix just means that the function input is expected to have a Point object as its first argument
    // since we didn't actually use p here, we could have written: ... = forall<double>([] PROTO_LAMBDA(Var<double, 1> v5, ...))


    // Here, the inputs D3 and D4 didn't have the same Box domain. Forall automatically computes the function on the union of the two domains:
    EXPECT_EQ(D5.box(), D3.box() & D4.box());

    //=============================================================================================
    // STENCIL

    // Unlike with forall, Stencils are standalone object functors. They are built using the Shift object which has the same API as Point
    Stencil<double> S = 0.5*Shift::X(-1) + 0.5*Shift::X(+1); // This means: S[D_{i}] = 0.5*D_{i-1} + 0.5*D_{i+1}
    BoxData<double, DIM> SD = S(D); // Apply the Stencil to create a new BoxData

    // You can also increment an existing BoxData with the output of the Stencil:
    SD += S(D);
    // Or replace the data in an existing BoxData with theoutput of the Stencil:
    SD |= S(D);

    // Stencil Objects conveniently know a lot of information about their input domain and output range. If you know one, you can compute the other:
    Box S0 = Box::Cube(32);
    Box SDomain = S.domain(S0); // The domain will be larger in the X direction since the Stencil needs ghost data along that coordinate
    Box SRange = S.range(S0); // Inversely, the range will be shrunk by the same amount

    EXPECT_EQ(SDomain, S0.grow(Point::X()));
    EXPECT_EQ(SRange, S0.grow(-Point::X()));

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
