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
