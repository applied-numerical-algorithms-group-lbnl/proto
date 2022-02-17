#include <gtest/gtest.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <cmath>
#include "Proto.H"
#include "UnitTestFunctions.H"

#define _USE_MATH_DEFINES

using namespace Proto;
using namespace std;
using namespace prototest;

int a_errorcode = 0;
bool a_didtestpass = true;

TEST(UnitTest, pointTest) {
    pointTest(a_errorcode, a_didtestpass);
    EXPECT_EQ(a_errorcode,0);
}

TEST(UnitTest, bxTest) {
    bxTest(a_errorcode, a_didtestpass);
    EXPECT_EQ(a_errorcode,0);
}

TEST(UnitTest, boxdataTest) {
    boxdataTest(a_errorcode, a_didtestpass);
    EXPECT_EQ(a_errorcode,0);
}

TEST(UnitTest, stencilTest) {
    stencilTest(a_errorcode, a_didtestpass);
    EXPECT_EQ(a_errorcode,0);
}

TEST(UnitTest, interpTest) {
    interpTest(a_errorcode, a_didtestpass);
    EXPECT_EQ(a_errorcode,0);
}

TEST(UnitTest, memtypeTest) {
    memtypeTest(a_errorcode, a_didtestpass);
    EXPECT_EQ(a_errorcode,0);
}
int main(int argc, char *argv[]) {
    int result = 0;
    ::testing::InitGoogleTest(&argc, argv);
#ifdef PR_MPI
    MPI_Init(&argc, &argv);
#endif
    result = RUN_ALL_TESTS();
#ifdef PR_MPI
    MPI_Finalize();
#endif
    return result;
}
