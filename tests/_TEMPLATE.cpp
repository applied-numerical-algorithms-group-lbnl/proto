#include <gtest/gtest.h>
#include "Proto.H"

using namespace Proto;

TEST(SUITE_NAME, TEST_NAME) {
    #if 0
    // use this to toggle tests on and off while developing.
    writeColorText(RED, "TEST DISABLED | This test should be either fixed or deleted");
    #else
    #if PR_VERBOSE > 0
    HDF5Handler h5;
    #endif
    // WRITE TEST CODE
    bool booleanTest = true;
    bool A = 1;
    bool B = 1;
    EXPECT_TRUE(booleanTest);
    EXPECT_EQ(A, B);
    #endif
}

// TEST...
// TEST...
// TEST...
// TEST...

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
