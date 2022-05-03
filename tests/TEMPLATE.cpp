#include <gtest/gtest.h>
#include "Proto.H"

using namespace Proto;

TEST(SUITE_NAME, TEST_NAME) {
    // WRITE TEST CODE
    bool booleanTest = true;
    bool A = 1;
    bool B = 1;
    EXPECT_TRUE(booleanTest);
    EXPECT_EQ(A, B);
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
