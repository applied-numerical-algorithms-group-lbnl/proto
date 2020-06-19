#include "gtest/gtest.h"
#include "Proto.H"
#include "Proto_Box.H"

using namespace Proto;

TEST(BoxTest, DefaultConstruction) {
  Box B;
  EXPECT_EQ(B.low(), Point::Zeros());
  EXPECT_EQ(B.high(), Point::Ones(-1));
  EXPECT_LE(B.size(), 0);
  EXPECT_TRUE(B.empty());
}

TEST(BoxTest, PointConstruction) {
  Point p1(1, 2, 3, 4, 5, 6);
  Point p0 = (-1) * p1;
  Box B0(p1, p0);
  Box B1(p0, p1);
  Box B2(p1);
  int s1 = 1;
  int s2 = 1;
  for (int ii = 0; ii < DIM; ii++) {
    s1 *= (p1[ii] - p0[ii] + 1);
    s2 *= (p1[ii] + 1);
  }
  EXPECT_LE(B0.size(), 0);
  EXPECT_TRUE(B0.empty());
  EXPECT_EQ(B1.size(), s1);
  EXPECT_FALSE(B1.empty());
  EXPECT_EQ(B2.size(), s2);
  EXPECT_FALSE(B2.empty());
}

TEST(BoxTest, CopyConstruction) {
  Box B0(Point(1, 2, 3, 4, 5, 6));
  Box B1(B0);
  EXPECT_EQ(B0, B1);
  EXPECT_NE(&B0,  &B1);
}

TEST(BoxTest, StaticMethods) {
  int size = 17;
  Box B = Box::Cube(size);
  EXPECT_EQ(B.size(), ipow<DIM>(size));
  EXPECT_EQ(B.low(), Point::Zeros());
  EXPECT_EQ(B.high(), Point::Ones(size - 1));
}

TEST(BoxTest, IterationAndIndexing) {
  Point max = Point(1, 2, 3, 4, 5, 6);
  Point min = -1 * max;
  Box B(min, max);
  int index = 0;
  Point p = min;
  for (auto iter = B.begin(); iter != B.end(); ++iter, ++index) {
    EXPECT_EQ(*iter, p);
    EXPECT_EQ(B.index(*iter), index);
    ++p[0];
    for (int ii = 0; ii < DIM - 1; ii++) {
      if (p[ii] > max[ii]) {
        p[ii] = min[ii];
        ++p[ii + 1];
      }
    }
  }
  p = max;
  index = B.size() - 1;
  for (auto iter = B.rbegin(); iter != B.rend(); --iter, --index) {
    EXPECT_EQ(*iter, p);
    EXPECT_EQ(B.index(*iter), index);
    --p[0];
    for (int ii = 0; ii < DIM - 1; ii++) {
      if (p[ii] < min[ii]) {
        p[ii] = max[ii];
        --p[ii + 1];
      }
    }
  }
}

TEST(BoxTest, Shift) {
  Box b0 = Box::Cube(16);
  Point s(1, -2, 3, -4, 5, -6);

  Box b1 = b0.shift(s);

  EXPECT_EQ(b0.low() + s, b1.low());
  EXPECT_EQ(b0.high() + s, b1.high());
  EXPECT_NE(&b0, &b1);
  EXPECT_EQ(b0, Box::Cube(16));
}

TEST(BoxTest, Grow) {
  Box b0 = Box::Cube(16);
  Point s(1, -2, 3, -4, 5, -6);
  //Grow (Point)
  Box b1 = b0.grow(s);

  EXPECT_EQ(b0.low() - s, b1.low());
  EXPECT_EQ(b0.high() + s, b1.high());
  EXPECT_NE(&b0, &b1);
  EXPECT_EQ(b0, Box::Cube(16));

  //Grow (scalar)
  b1 = b0.grow(3);

  EXPECT_EQ(b0.low() - 3, b1.low());
  EXPECT_EQ(b0.high() + 3, b1.high());
  EXPECT_NE(&b0, &b1);
  EXPECT_EQ(b0, Box::Cube(16));
}

TEST(BoxTest, Coarsen) {
  Point low = Point::Ones(-2);
  Point high = Point::Ones(3);
  Point r = Point::Ones() + Point::Basis(0); //(2,1,1,...,1)
  Box b0 = Box(low, high);
  Box b1 = b0.coarsen(2);
  Box b2 = b0.coarsen(3);
  Box b3 = b0.coarsen(r);

  EXPECT_NE(&b0, &b1);
  EXPECT_NE(&b0, &b2);
  EXPECT_NE(&b0, &b3);

  EXPECT_EQ(b1.low(), b0.low() / 2);
  EXPECT_EQ(b1.high(), b0.high() / 2);
  EXPECT_EQ(b2.low(), b0.low() / 3);
  EXPECT_EQ(b2.high(), b0.high() / 3);
  for (int ii = 0; ii < DIM; ++ii) {
    EXPECT_EQ(b3.low()[ii], b0.low()[ii]/r[ii]);
    EXPECT_EQ(b3.high()[ii], b0.high()[ii]/r[ii]);
  }
}

TEST(BoxTest, Refine) {
  Point low = Point::Ones(-2);
  Point high = Point::Ones(3);
  Point r = Point::Ones() + Point::Basis(0); //(2,1,1,...,1)
  Box b0 = Box(low, high);
  Box b1 = b0.refine(2);
  Box b2 = b0.refine(3);
  Box b3 = b0.refine(r);

  EXPECT_NE(&b0, &b1);
  EXPECT_NE(&b0, &b2);
  EXPECT_NE(&b0, &b3);

  EXPECT_EQ(b1.low(), b0.low() * 2);
  EXPECT_EQ(b1.high(), (b0.high() + Point::Ones()) * 2 - Point::Ones());
  EXPECT_EQ(b2.low(), b0.low() * 3);
  EXPECT_EQ(b2.high(), (b0.high() + Point::Ones()) * 3 - Point::Ones());
  EXPECT_EQ(b3.low()[0], b0.low()[0] * 2);
  EXPECT_EQ(b3.high()[0], (b0.high()[0] + 1) * 2 - 1);

  for (int ii = 1; ii < DIM; ii++) {
    EXPECT_EQ(b3.low()[ii], b0.low()[ii]);
    EXPECT_EQ(b3.high()[ii], b0.high()[ii]);
  }
}

TEST(BoxTest, Flatten) {
  Box b0 = Box::Cube(17);
  for (int ii = 0; ii < DIM; ii++) {
    Box b1 = b0.flatten(ii);
    Box b2 = b0.flatten(ii, true);
    EXPECT_EQ(b1.low(), b0.low());
    EXPECT_EQ(b2.high(), b0.high());
    for (int jj = 0; jj < DIM; jj++) {
      if (jj == ii) {
        EXPECT_EQ(b1.high()[jj], b1.low()[jj]);
        EXPECT_EQ(b2.low()[jj], b2.high()[jj]);
      } else {
        EXPECT_EQ(b1.high()[jj], b0.high()[jj]);
        EXPECT_EQ(b2.low()[jj], b0.low()[jj]);
      }
    }
  }
}

TEST(BoxTest, Extrude) {
  Box b0 = Box::Cube(17);
  for (int ii = 0; ii < DIM; ii++) {
    Box b1 = b0.extrude(ii, 3, false);
    Box b2 = b0.extrude(ii, 3, true);
    Box b3 = b0.extrude(ii, 3);
    EXPECT_EQ(b1.high(), b0.high());
    EXPECT_EQ(b2.low(), b0.low());
    for (int jj = 0; jj < DIM; jj++) {
      if (jj == ii) {
        EXPECT_EQ(b1.low()[jj], b0.low()[jj] - 3);
        EXPECT_EQ(b2.high()[jj], b0.high()[jj] + 3);
      } else {
        EXPECT_EQ(b1.low()[jj], b0.low()[jj]);
        EXPECT_EQ(b2.high()[jj], b0.high()[jj]);
      }
    }
    EXPECT_EQ(b3, b2);
  }
}

TEST(BoxTest, Mod) {
  Point high(1, 2, 3, 4, 5, 6);
  Point low = high * (-1);
  Box B(low, high);

  Point p0 = Point::Zeros();
  Point p1(1, 2, 3, 4, 5, 6);
  Point p2 = high + Point::Ones();
  Point p3 = low - Point::Ones();
  Point p4 = B.flatten(0).high() - Point::Basis(0);
  Point p5 = B.flatten(0, true).low() + Point::Basis(0);

  Point q0 = B.mod(p0);
  Point q1 = B.mod(p1);
  Point q2 = B.mod(p2);
  Point q3 = B.mod(p3);
  Point q4 = B.mod(p4);
  Point q5 = B.mod(p5);

  EXPECT_EQ(q0, p0);
  EXPECT_EQ(q1, p1);
  EXPECT_EQ(q2, low);
  EXPECT_EQ(q3, high);
  EXPECT_EQ(q4, high);
  EXPECT_EQ(q5, low);
}
