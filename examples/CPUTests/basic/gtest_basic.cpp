#include "gtest/gtest.h"

TEST (gtest_basic, test_eq) {
  EXPECT_EQ(1, 1);
}

TEST (gtest_basic, test_neq) {
  EXPECT_NE(1, 0);
}
