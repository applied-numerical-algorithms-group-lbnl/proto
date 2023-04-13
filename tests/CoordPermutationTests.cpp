#include <gtest/gtest.h>
#include "Proto.H"
#include "Lambdas.H"

using namespace Proto;

void compareMatrices(
        Array<Array<int, DIM>, DIM>& M0,
        Array<Array<int, DIM>, DIM>& M1)
{
    for (int ii = 0; ii < DIM; ii++)
    {
        for (int jj = 0; jj < DIM; jj++)
        {
            EXPECT_EQ(M0[ii][jj], M1[ii][jj]);
        }
    }
}

TEST(CoordPermutation, Identity)
{
    auto I = CoordPermutation::identity();
    auto M = I.matrix();
#if DIM==2
    Array<Array<int, DIM>, DIM> M0{{{{1,0}},{{0,1}}}};
#elif DIM==3
    Array<Array<int, DIM>, DIM> M0{{{{1,0,0}},{{0,1,0}},{{0,0,1}}}};
#else
    MayDay<void>::Abort("CoordPermutation test not written for this value of DIM");
#endif
    compareMatrices(M, M0);
}

TEST(CoordPermutation, Clockwise)
{
    auto R = CoordPermutation::cw();
    auto M = R.matrix();
    auto basis = Point::Ones().parallelUnit();
#if DIM==2
    Array<Array<int, DIM>, DIM> M0{{{{0,-1}},{{1,0}}}};
#elif DIM==3
    Array<Array<int, DIM>, DIM> M0{{{{0,-1,0}},{{1,0,0}},{{0,0,1}}}};
    for (int ii = 0; ii < 3; ii++)
    {
        auto Ri = CoordPermutation::cw(ii);
        Point p0 = basis[(ii+1)%3];
        Point p1 = basis[(ii+2)%3];
        Point q0 = Ri(p0);
        Point q1 = Ri(p1);
        EXPECT_EQ(q0, p1);
        EXPECT_EQ(q1, -p0);
    }
#else
    MayDay<void>::Abort("CoordPermutation test not written for this value of DIM");
#endif
    compareMatrices(M, M0);
}

TEST(CoordPermutation, CounterClockwise)
{
    auto R = CoordPermutation::ccw();
    auto M = R.matrix();
    auto basis = Point::Ones().parallelUnit();
#if DIM==2
    Array<Array<int, DIM>, DIM> M0{{{{0,1}},{{-1,0}}}};
#elif DIM==3
    Array<Array<int, DIM>, DIM> M0{{{{0,1,0}},{{-1,0,0}},{{0,0,1}}}};
    for (int ii = 0; ii < 3; ii++)
    {
        auto Ri = CoordPermutation::ccw(ii);
        Point p0 = basis[(ii+1)%3];
        Point p1 = basis[(ii+2)%3];
        Point q0 = Ri(p0);
        Point q1 = Ri(p1);
        EXPECT_EQ(q0, -p1);
        EXPECT_EQ(q1, p0);
    }
#else
    MayDay<void>::Abort("CoordPermutation test not written for this value of DIM");
#endif
    compareMatrices(M, M0);
}

TEST(CoordPermutation, Reverse)
{
    auto R = CoordPermutation::reverse();
    auto M = R.matrix();
#if DIM==2
    Array<Array<int, DIM>, DIM> M0{{{{-1,0}},{{0,-1}}}};
#elif DIM==3
    Array<Array<int, DIM>, DIM> M0{{{{-1,0,0}},{{0,-1,0}},{{0,0,-1}}}};
#else
    MayDay<void>::Abort("CoordPermutation test not written for this value of DIM");
#endif
    compareMatrices(M, M0);
}

TEST(CoordPermutation, InvertAxis)
{
    CoordPermutation R({{1,1,-1}});
    auto M = R.matrix();
#if DIM==2
    Array<Array<int, DIM>, DIM> M0{{{{1,0}},{{0,-1}}}};
#elif DIM==3
    Array<Array<int, DIM>, DIM> M0{{{{1,0,0}},{{0,-1,0}},{{0,0,1}}}};
#else
    MayDay<void>::Abort("CoordPermutation test not written for this value of DIM");
#endif
    compareMatrices(M, M0);
}

TEST(CoordPermutation, RotatePlane)
{
    auto R01 = CoordPermutation::rotatePlane(0,1);
    auto R10 = CoordPermutation::rotatePlane(1,0);
    auto CW = CoordPermutation::cw();
    auto CCW = CoordPermutation::ccw();
    EXPECT_EQ(R01, CW);
    EXPECT_EQ(R10, CCW);
}

TEST(CoordPermutation, Apply2D)
{
    Point p(7,-6,5,-4,3,-2,1);
    auto I = CoordPermutation::identity();
    auto R = CoordPermutation::reverse();
    auto CW = CoordPermutation::cw();
    auto CCW = CoordPermutation::ccw();
    Point p_I = I(p);
    Point p_R = R(p);
    Point p_CW = CW(p);
    Point p_CCW = CCW(p);
    for (int ii = 0; ii < DIM; ii++)
    {
        EXPECT_EQ(p[ii], p_I[ii]);
    }
    EXPECT_EQ(p[0], -p_R[0]);
    EXPECT_EQ(p[1], -p_R[1]);
    EXPECT_EQ(p[0], p_CW[1]);
    EXPECT_EQ(p[1], -p_CW[0]);
    EXPECT_EQ(p[0], -p_CCW[1]);
    EXPECT_EQ(p[1], p_CCW[0]);
#if DIM > 2
    EXPECT_EQ(p[2], -p_R[2]);
    EXPECT_EQ(p[2], p_CW[2]);
    EXPECT_EQ(p[2], p_CCW[2]);
#endif
}

TEST(CoordPermutation, Variadic)
{
    Array<int, DIM> mapping;
    Array<int, DIM> mirror;
    mirror.fill(1);
    for (int ii = 0; ii < DIM; ii++)
    {
        mapping[ii] = (ii+1)%DIM;
        if (ii % 2 == 1) { mirror[ii] = -1; }
    }
#if DIM == 2
    CoordPermutation R{
        {0, mapping[0], mirror[0]},
        {1, mapping[1], mirror[1]}};
#elif DIM == 3
    CoordPermutation R{
        {0, mapping[0], mirror[0]},
        {1, mapping[1], mirror[1]},
        {2, mapping[2], mirror[2]}};
#else
    MayDay<void>::Abort("CoordPermutation test not written for this value of DIM");
#endif
    Point p0(0,-1,2,-3,4,-5,6);
    Point p1 = R(p0);
    EXPECT_EQ(p0[0], p1[mapping[0]]*mirror[0]);
    EXPECT_EQ(p0[1], p1[mapping[1]]*mirror[1]);
#if DIM > 2
    EXPECT_EQ(p0[2], p1[mapping[2]]*mirror[2]);
#endif
}

TEST(CoordPermutation, VariadicPoint)
{
    Point x = Point::Basis(0);
    Point y = Point::Basis(1);
    CoordPermutation R{{x,y},{y,-x}};
    auto Rx = R(x);
    auto Ry = R(y);
    EXPECT_EQ(Rx,  y);
    EXPECT_EQ(Ry, -x);
}

TEST(CoordPermutation, Inverse)
{
    auto CW = CoordPermutation::cw();
    auto CCW = CoordPermutation::ccw();
    EXPECT_EQ(CW, CCW.inverse());
    EXPECT_EQ(CCW, CW.inverse());
}

TEST(CoordPermutation, Convolution)
{
    Point p0(0,-1,2,-3,4,-5,6);
    auto CW = CoordPermutation::cw();
    auto CCW = CoordPermutation::ccw();

    auto CW2 = CW*CW;
    auto CW3 = CW2*CW;
    auto CCW2 = CCW*CCW;
    auto CCW3 = CCW2*CCW;
    
    EXPECT_EQ(CW, CCW3);
    EXPECT_EQ(CW2, CCW2);
    EXPECT_EQ(CW3, CCW);

    auto pCW = CW(p0);
    auto pCW2 = CW2(p0);
    auto pCW3 = CW3(p0);
    auto pCCW = CCW(p0);
    auto pCCW2 = CCW2(p0);
    auto pCCW3 = CCW3(p0);


    EXPECT_EQ(pCW, pCCW3);
    EXPECT_EQ(pCW2, pCCW2);
    EXPECT_EQ(pCW3, pCCW);
}

TEST(CoordPermutation, RotatePoint)
{
    auto CW = CoordPermutation::cw();
    auto R = CoordPermutation{{0,0,-1}};
    int dx = 4;
    int dy = 5;
    Box B(Point(dx,dy,6,7,8,9));
    B = B.shift(Point::Ones());
    Box B_CW(Point(dy,dx,6,7,8,9));
    B_CW = B_CW.shift(Point::Ones());
    Box B_R(Point(dx,dy,6,7,8,9));
    B_R = B_R.shift(Point::Ones());
    Point s_cw = Point::Basis(0, dy) + B_CW.low();
    Point s_r = Point::Basis(0, dx) + B_R.low();
    for (auto pi : B)
    {

        auto p_cw = CW.rotatePoint(pi, B, B_CW);
        auto p_r = R.rotatePoint(pi, B, B_R);
        EXPECT_EQ(p_cw, CW(pi - B.low()) + s_cw);
        EXPECT_EQ(p_r, R(pi - B.low()) + s_r);
    }
}

TEST(CoordPermutation, RotateCell)
{
    auto CW = CoordPermutation::cw();
    auto R = CoordPermutation{{0,0,-1}};
    int dx = 4;
    int dy = 5;
    Box B(Point(dx,dy,6,7,8,9));
    B = B.shift(Point::Ones());
    Box B_CW(Point(dy,dx,6,7,8,9));
    B_CW = B_CW.shift(Point::Ones());
    Box B_R(Point(dx,dy,6,7,8,9));
    B_R = B_R.shift(Point::Ones());
    Point s_cw = Point::Basis(0, dy-1) + B_CW.low();
    Point s_r = Point::Basis(0, dx-1) + B_R.low();
    for (auto pi : B)
    {
        auto p_cw = CW.rotateCell(pi, B, B_CW);
        auto p_r = R.rotateCell(pi, B, B_R);
        EXPECT_TRUE(B_CW.contains(p_cw));
        EXPECT_TRUE(B_R.contains(p_r));
        EXPECT_EQ(p_cw, CW(pi - B.low()) + s_cw);
        EXPECT_EQ(p_r, R(pi - B.low()) + s_r);
    }
}
TEST(CoordPermutation, RotateBuffer)
{
    constexpr unsigned int C = 3;
    auto CW = CoordPermutation::cw();
    auto R = CoordPermutation{{0,0,-1}};
    int dx = 4;
    int dy = 5;
    Box B0(Point(dx,dy,6,7,8,9));
    B0 = B0.shift(Point::Ones());
    Box B_CW(Point(dy,dx,6,7,8,9));
    B_CW = B_CW.shift(Point::Ones());
    Box B_R(Point(dx,dy,6,7,8,9));
    B_R = B_R.shift(Point::Ones());
    BoxData<int, C, HOST> hostSrc(B0);
    BoxData<int, C, HOST> hostDst_CW(B0);
    BoxData<int, C, HOST> hostDst_R(B0);
    forallInPlace_p(f_pointID, hostSrc);
    forallInPlace_p(f_pointID, hostDst_CW);
    forallInPlace_p(f_pointID, hostDst_R);
    hostDst_CW.rotate(B_CW, CW);
    hostDst_R.rotate(B_R, R);
    
    EXPECT_EQ(hostDst_CW.box(), B_CW);
    EXPECT_EQ(hostDst_R.box(), B_R);
    for (int cc = 0; cc < C; cc++)
    {
        for (auto pi : B0)
        {
            auto p_CW = CW.rotateCell(pi, B0, B_CW);
            auto p_R  = R.rotateCell(pi, B0, B_R);
            EXPECT_EQ(hostSrc(pi, cc), hostDst_R(p_R, cc));
            EXPECT_EQ(hostSrc(pi, cc), hostDst_CW(p_CW, cc));
        }
    }

    hostDst_CW.define(B_CW);
    hostDst_CW.setVal(0);
    hostDst_R.define(B_R);
    hostDst_R.setVal(0);
    hostSrc.copyTo(hostDst_CW, CW);
    hostSrc.copyTo(hostDst_R, R);
    
    for (int cc = 0; cc < C; cc++)
    {
        for (auto pi : B0)
        {
            auto p_CW = CW.rotateCell(pi, B0, B_CW);
            auto p_R  = R.rotateCell(pi, B0, B_R);
            EXPECT_EQ(hostSrc(pi, cc), hostDst_R(p_R, cc));
            EXPECT_EQ(hostSrc(pi, cc), hostDst_CW(p_CW, cc));
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
