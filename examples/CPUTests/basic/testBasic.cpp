#include "gtest/gtest.h"

TEST (TestBasic, TestEq) {
  EXPECT_EQ(1, 1) << "One should always equal to one";
  EXPECT_TRUE(true) << "True is right";
  EXPECT_FALSE(false) << "False is wrong";
}

TEST (TestBasic, TestNEq) {
  EXPECT_NE(1, 0) << "One should not equal to zero";
  EXPECT_NE(true, false) << "True should not be false";
}

TEST (TestBasic, TwosComplement) {
  EXPECT_EQ((long)(-1), (long)(~0ul)) << "True should not be false";
}

TEST (TestBasic, IntegerModWMask) {
  int MOD = 16;
  int mask = MOD - 1;
  int mask_r = -MOD;
  EXPECT_EQ(((6 * MOD + 5) & mask), 5);
  EXPECT_EQ(((-6 * MOD + 5) & mask), 5);
  EXPECT_EQ(((6 * MOD + 5) & mask_r), (6 * MOD));
  EXPECT_EQ(((-6 * MOD + 5) & mask_r), (-6 * MOD));
}