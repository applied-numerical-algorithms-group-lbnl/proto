#include <gtest/gtest.h>
#include "ProtoMMB.H"

using namespace Proto;

void compareMatrices(
        std::array<std::array<int, DIM>, DIM>& M0,
        std::array<std::array<int, DIM>, DIM>& M1)
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
    std::array<std::array<int, DIM>, DIM> M0{{{{1,0}},{{0,1}}}};
#elif DIM==3
    std::array<std::array<int, DIM>, DIM> M0{{{{1,0,0}},{{0,1,0}},{{0,0,1}}}};
#else
    MayDay<void>::Abort("CoordPermutation test not written for this value of DIM");
#endif
    compareMatrices(M, M0);
}

TEST(CoordPermutation, Clockwise)
{
    auto R = CoordPermutation::cw();
    auto M = R.matrix();
#if DIM==2
    std::array<std::array<int, DIM>, DIM> M0{{{{0,-1}},{{1,0}}}};
#elif DIM==3
    std::array<std::array<int, DIM>, DIM> M0{{{{0,-1,0}},{{1,0,0}},{{0,0,1}}}};
#else
    MayDay<void>::Abort("CoordPermutation test not written for this value of DIM");
#endif
    compareMatrices(M, M0);
}

TEST(CoordPermutation, CounterClockwise)
{
    auto R = CoordPermutation::ccw();
    auto M = R.matrix();
#if DIM==2
    std::array<std::array<int, DIM>, DIM> M0{{{{0,1}},{{-1,0}}}};
#elif DIM==3
    std::array<std::array<int, DIM>, DIM> M0{{{{0,1,0}},{{-1,0,0}},{{0,0,1}}}};
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
    std::array<std::array<int, DIM>, DIM> M0{{{{-1,0}},{{0,-1}}}};
#elif DIM==3
    std::array<std::array<int, DIM>, DIM> M0{{{{-1,0,0}},{{0,-1,0}},{{0,0,-1}}}};
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
    std::array<int, DIM> mapping;
    std::array<int, DIM> mirror;
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
    int dx = 4;
    int dy = 5;
    Box B(Point(dx,dy,6,7,8,9));
    auto CW = CoordPermutation::cw();
    auto R = CoordPermutation{{0,0,-1}};
    Point s_cw = Point::Basis(0, dy);
    Point s_r = Point::Basis(0, dx);
    for (auto pi : B)
    {
        auto p_cw = CW.rotatePoint(pi, B);
        auto p_r = R.rotatePoint(pi, B);
        EXPECT_EQ(p_cw, CW(pi) + s_cw);
        EXPECT_EQ(p_r, R(pi) + s_r);
    }
}

TEST(CoordPermutation, RotateCell)
{
    int dx = 4;
    int dy = 5;
    Box B(Point(dx,dy,6,7,8,9));
    auto CW = CoordPermutation::cw();
    auto R = CoordPermutation{{0,0,-1}};
    Point s_cw = Point::Basis(0, dy-1);
    Point s_r = Point::Basis(0, dx-1);
    for (auto pi : B)
    {
        auto p_cw = CW.rotateCell(pi, B);
        auto p_r = R.rotateCell(pi, B);
        EXPECT_EQ(p_cw, CW(pi) + s_cw);
        EXPECT_EQ(p_r, R(pi) + s_r);
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
