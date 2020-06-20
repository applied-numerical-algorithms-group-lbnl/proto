#include "gtest/gtest.h"
#include "Proto.H"
#include "Proto_BoxData.H"

using namespace Proto;

#ifdef PROTO_BRICK

#include "brick.h"
#include "bricksetup.h"
#include "multiarray.h"

namespace Proto {
  BrickMetaCollector brickMetaCollector;
}

#define ValType bElem
#else
#define ValType int
#endif

TEST(BoxDataTest, DefaultConstruction) {
  BoxData<ValType, 2, 3> BD;

  EXPECT_EQ(BD.box(), Box(Point::Ones(-1)));
  EXPECT_EQ(BD.box(), Box());
  EXPECT_EQ(BD.padded_box(), BD.box());
  EXPECT_LE(BD.size(), 0);
  EXPECT_LE(BD.padded_size(), 0);
}

TEST(BoxDataTest, BoxConstructor) {
  Point p(31, 63, 31, 63, 127);
  Box B = Box(p);
  BoxData<ValType, 3, 4, 5> BD(B);
  int size = 1;
  for (int ii = 0; ii < DIM; ii++)
    size *= p[ii] + 1;
  EXPECT_EQ(BD.box(), B);
  EXPECT_EQ(BD.size(), size * 3 * 4 * 5);
}

TEST(BoxDataTest, BoxConstructorStack) {
  Point p(31, 63, 31, 63, 127);
  Box B = Box(p);
  BoxData<ValType, 3, 4, 5> BD(B, true);
  int size = 1;
  for (int ii = 0; ii < DIM; ii++)
    size *= p[ii] + 1;
  EXPECT_EQ(BD.box(), B);
  EXPECT_EQ(BD.size(), size * 3 * 4 * 5);
}

#ifdef PROTO_BRICK
TEST(BoxDataTest, BoxConstructorBrickFail) {
  EXPECT_THROW(({
    BoxData<ValType, 3>(Box(Point(32, 32, 32, 32, 32)));
  }), std::runtime_error) << "Brick cannot construct Boxes with odd dimensions (33, ...)";
}

#endif

TEST(BoxDataTest, InitConstructorWithAccessor) {
  Point p(31, 63, 31, 63, 127);
  Box B = Box(p);
  ValType val = static_cast<ValType>(1337);
  BoxData<ValType, 3, 4, 5> BD(B, val);
  int size = 1;
  for (int ii = 0; ii < DIM; ii++)
    size *= p[ii] + 1;
  EXPECT_EQ(BD.box(), B);
  EXPECT_EQ(BD.size(), size * 3 * 4 * 5);
  ValType maxval = BD.max();
  ValType minval = BD.min();
  EXPECT_EQ(maxval, val);
  EXPECT_EQ(minval, val);
}

class ForallTest : public ::testing::Test {
protected:
  ForallTest() : p(31, 63, 31, 63, 127), box(p), boxData(box, (ValType) 0) {}

  Point p;
  Box box;
  BoxData<ValType, DIM> boxData;
};

TEST_F(ForallTest, forall) {
  auto setMagic = [](Var<ValType, DIM> &var) {
    for (int i = 0; i < DIM; ++i)
      var(i) = i + 1;
  };
  // Partial set
  Point p_partial = Point(15, 31, 15, 31, 63);
  forallInPlace(setMagic, Box(p_partial), boxData);
  {
    auto first = boxData.var(Point::Zeros());
    auto mid = boxData.var(p_partial);
    auto midp = boxData.var(p_partial + Point::Ones());
    auto last = boxData.var(p);
    for (int i = 0; i < DIM; ++i) {
      EXPECT_EQ(first(i), i + 1);
      EXPECT_EQ(mid(i), first(i));
      EXPECT_EQ(midp(i), 0);
      EXPECT_EQ(last(i), 0);
    }
  }
  // Partial set out-of-place
  // This uses stack allocator
  auto mulMagic = [](Var<ValType, DIM> &dst, Var<ValType, DIM> &src) {
    for (int i = 0; i < DIM; ++i)
      dst(i) = (i + 1) * src(i);
  };
  auto forall_dst = forall<ValType, DIM>(mulMagic, Box(p_partial), boxData);
  {
    EXPECT_EQ(forall_dst.box(), Box(p_partial));
    EXPECT_EQ(forall_dst.box(), forall_dst.padded_box());
    auto first = forall_dst.var(Point::Zeros());
    auto last = forall_dst.var(p_partial);
    auto orig_first = boxData.var(Point::Zeros());
    auto orig_last = boxData.var(p);
    for (int i = 0; i < DIM; ++i) {
      EXPECT_EQ(first(i), (i + 1) * (i + 1));
      EXPECT_EQ(last(i), first(i));
      EXPECT_EQ(orig_first(i), i + 1);
      EXPECT_EQ(orig_last(i), 0);
    }
  }
  // Intersect
  Box intersect = Box(p_partial).shift(Point::Ones());
  forallInPlace(setMagic, intersect, forall_dst);
  {
    auto first = forall_dst.var(Point::Zeros());
    auto mid = forall_dst.var(Point::Ones());
    auto last = forall_dst.var(p_partial);
    for (int i = 0; i < DIM; ++i) {
      EXPECT_EQ(first(i), (i + 1) * (i + 1));
      EXPECT_EQ(mid(i), i + 1);
      EXPECT_EQ(last(i), mid(i));
    }
  }
  // Empty box
  forallInPlace(setMagic, Box(), forall_dst);
  {
    auto first = forall_dst.var(Point::Zeros());
    auto last = forall_dst.var(p_partial);
    for (int i = 0; i < DIM; ++i) {
      EXPECT_EQ(first(i), i + 1);
      EXPECT_EQ(last(i), first(i));
    }
  }
}

TEST_F(ForallTest, forall_p) {
  auto setMagic = [](const Point &p, Var<ValType, DIM> &var) {
    for (int i = 0; i < DIM; ++i)
      var(i) = p[i];
  };
  // Partial set
  Point p_partial = Point(15, 31, 15, 31, 63);
  Point p_arb = Point(1, 2, 3, 4, 5, 6);
  forallInPlace_p(setMagic, Box(p_partial), boxData);
  {
    auto first = boxData.var(p_arb);
    auto mid = boxData.var(p_partial);
    auto midp = boxData.var(p_partial + Point::Ones());
    auto last = boxData.var(p);
    for (int i = 0; i < DIM; ++i) {
      EXPECT_EQ(first(i), i + 1);
      EXPECT_EQ(mid(i), p_partial[i]);
      EXPECT_EQ(midp(i), 0);
      EXPECT_EQ(last(i), 0);
    }
  }
  // Partial set out-of-place
  // This uses stack allocator
  auto mulMagic = [](const Point &p, Var<ValType, DIM> &dst, Var<ValType, DIM> &src) {
    for (int i = 0; i < DIM; ++i)
      dst(i) = p[i] * src(i);
  };
  auto forall_dst = forall_p<ValType, DIM>(mulMagic, Box(p_partial), boxData);
  {
    EXPECT_EQ(forall_dst.box(), Box(p_partial));
    EXPECT_EQ(forall_dst.box(), forall_dst.padded_box());
    auto first = forall_dst.var(p_arb);
    auto last = forall_dst.var(p_partial);
    auto orig_first = boxData.var(p_arb);
    auto orig_last = boxData.var(p);
    for (int i = 0; i < DIM; ++i) {
      EXPECT_EQ(first(i), (i + 1) * (i + 1));
      EXPECT_EQ(last(i), p_partial[i] * p_partial[i]);
      EXPECT_EQ(orig_first(i), i + 1);
      EXPECT_EQ(orig_last(i), 0);
    }
  }
  // Intersect
  Box intersect = Box(p_partial).shift(Point::Ones());
  forallInPlace_p(setMagic, intersect, forall_dst);
  {
    auto first = forall_dst.var(p_arb - Point::Ones());
    auto mid = forall_dst.var(p_arb + Point::Ones());
    auto last = forall_dst.var(p_partial);
    for (int i = 0; i < DIM; ++i) {
      EXPECT_EQ(first(i), i * i);
      EXPECT_EQ(mid(i), i + 2);
      EXPECT_EQ(last(i), p_partial[i]);
    }
  }
}