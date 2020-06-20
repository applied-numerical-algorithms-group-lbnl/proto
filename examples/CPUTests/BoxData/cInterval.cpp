#include "gtest/gtest.h"
#include "Proto.H"
#include "Proto_BoxData.H"

using namespace Proto;

TEST(CIntervalTest, Construction) {

  CInterval I0(1, 2, 3, 4, 5, 6);
  CInterval I1{{1, 2},
               {3, 4},
               {5, 6}};
  CInterval I2{{},
               {3, 4},
               {}};
  CInterval I3{1, 2};

  EXPECT_EQ(I0.low(0), 1);
  EXPECT_EQ(I0.high(0), 2);
  EXPECT_EQ(I0.low(1), 3);
  EXPECT_EQ(I0.high(1), 4);
  EXPECT_EQ(I0.low(2), 5);
  EXPECT_EQ(I0.high(2), 6);

  EXPECT_EQ(I1.low(0), 1);
  EXPECT_EQ(I1.high(0), 2);
  EXPECT_EQ(I1.low(1), 3);
  EXPECT_EQ(I1.high(1), 4);
  EXPECT_EQ(I1.low(2), 5);
  EXPECT_EQ(I1.high(2), 6);

  EXPECT_EQ(I2.low(0), 0);
  EXPECT_EQ(I2.high(0), 0);
  EXPECT_EQ(I2.low(1), 3);
  EXPECT_EQ(I2.high(1), 4);
  EXPECT_EQ(I2.low(2), 0);
  EXPECT_EQ(I2.high(2), 0);

  EXPECT_EQ(I3.low(0), 1);
  EXPECT_EQ(I3.high(0), 2);
  EXPECT_EQ(I3.low(1), 0);
  EXPECT_EQ(I3.high(1), 0);
  EXPECT_EQ(I3.low(2), 0);
  EXPECT_EQ(I3.high(2), 0);
}
