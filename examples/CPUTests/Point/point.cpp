#include "gtest/gtest.h"
#include "Proto_Point.H"

using namespace Proto;

TEST(PointTest, DefaultConstruction) {
  Point p;
  for (int ii = 0; ii < DIM; ii++)
    EXPECT_EQ(p[ii], 0);
}

TEST(PointTest, ArrayConstruction) {
  int v[DIM];
  for (int ii = 0; ii < DIM; ii++)
    v[ii] = (ii + 1) * 17;
  Point p(v);
  for (int ii = 0; ii < DIM; ii++)
    EXPECT_EQ(p[ii], v[ii]);
}

TEST(PointTest, VariadicConstruction) {
  Point p(1, 2, 3, 4, 5, 6);
  for (int ii = 0; ii < DIM; ii++)
    EXPECT_EQ(p[ii], ii + 1);
}

TEST(PointTest, CopyConstruction) {
  Point p(1, 2, 3, 4, 5, 6);
  Point q(p);
  for (int ii = 0; ii < DIM; ii++)
    EXPECT_EQ(p[ii], q[ii]);
}

TEST(PointTest, StaticMethods) {
  Point p0 = Point::Zeros();
  Point p1 = Point::Ones();
  Point p2 = Point::Ones(17);
  Point p3 = Point::Basis(0);
  Point p4 = Point::Basis(DIM - 1, 17);
  for (int ii = 0; ii < DIM; ii++) {
    EXPECT_EQ(p0[ii], 0);
    EXPECT_EQ(p1[ii], 1);
    EXPECT_EQ(p2[ii], 17);
    EXPECT_EQ(p3[ii], ii == 0 ? 1 : 0);
    EXPECT_EQ(p4[ii], ii == DIM - 1 ? 17 : 0);
  }
}

TEST(PointTest, AccessorMethods) {
  const Point q(1, 2, 3, 4, 5, 6);
  for (int ii = 0; ii < DIM; ii++)
    EXPECT_EQ(q[ii], q.m_tuple[ii]);
}

TEST(PointTest, AlgebricOperators) {
  Point p0(1, 2, 3, 4, 5, 6);
  Point p1;
  int flip = -1;
  for (int ii = 0; ii < DIM; ii++) {
    p1[ii] = p0[ii] * 17 * flip;
    flip *= -1;
  }

  Point p2 = p0 + p1;
  for (int ii = 0; ii < DIM; ii++)
    EXPECT_EQ(p2[ii], p0[ii] + p1[ii]);

  p2 = p1 - p0;
  for (int ii = 0; ii < DIM; ii++)
    EXPECT_EQ(p2[ii], p1[ii] - p0[ii]);

  p2 = p1 * p0;
  for (int ii = 0; ii < DIM; ii++)
    EXPECT_EQ(p2[ii], p1[ii] * p0[ii]);

  p2 = p1 / p0;
  for (int ii = 0; ii < DIM; ii++)
    EXPECT_EQ(p2[ii], p1[ii] / p0[ii]);

  p2 = p1 / 17;
  for (int ii = 0; ii < DIM; ii++)
    EXPECT_EQ(p2[ii], p1[ii] / 17);

  p1 = Point::Ones(2);
  p2 = p0 % p1;
  for (int ii = 0; ii < DIM; ii++)
    EXPECT_EQ(p2[ii], p0[ii] % p1[ii]);
}
