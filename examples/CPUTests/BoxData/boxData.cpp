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

TEST(BoxTest, DefaultConstruction) {
  BoxData<ValType, 2, 3> BD;

  EXPECT_EQ(BD.box(), Box(Point::Ones(-1)));
  EXPECT_EQ(BD.box(), Box());
  EXPECT_EQ(BD.padded_box(), BD.box());
  EXPECT_LE(BD.size(), 0);
  EXPECT_LE(BD.padded_size(), 0);
}

TEST(BoxTest, BoxConstructor) {
  Point p(31,63,31,63,127);
  Box B = Box(p);
  BoxData<ValType, 3, 4, 5> BD(B);
  int size = 1;
  for (int ii = 0; ii < DIM; ii++)
    size *= p[ii] + 1;
  EXPECT_EQ(BD.box(), B);
  EXPECT_EQ(BD.size(), size*3*4*5);
}

#ifdef PROTO_BRICK
TEST(BoxTest, BoxConstructorBrickFail) {
  EXPECT_THROW(({
    BoxData<ValType, 3>(Box(Point(32, 32, 32, 32, 32)));
  }), std::runtime_error) << "Brick cannot construct Boxes with odd dimensions (33, ...)";
}
#endif

TEST(BoxTest, InitConstructorWithAccessor) {
  Point p(31,63,31,63,127);
  Box B = Box(p);
  ValType val = static_cast<ValType>(1337);
  BoxData<ValType, 3, 4, 5> BD(B, val);
  int size = 1;
  for (int ii = 0; ii < DIM; ii++)
    size *= p[ii] + 1;
  EXPECT_EQ(BD.box(), B);
  EXPECT_EQ(BD.size(), size*3*4*5);
  ValType maxval = BD.max();
  ValType minval = BD.min();
  EXPECT_EQ(maxval, val);
  EXPECT_EQ(minval, val);
}
