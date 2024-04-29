#include <gtest/gtest.h>
#include "Proto.H"

using namespace Proto;

TEST(EXAMPLE, TEST_1) {
    // ...
    // normal development code
    // ...
    EXPECT_TRUE(boolean_test)
    EXPECT_EQ(a,b)
    EXPECT_LT(value, tolerance)
    // ...
    // normal development code
    // ...

#if PR_VERBOSE > 0

    // HDF5Handler h5;
    // h5.writeLevel(..., diagnosticData, ...);
    // pr_out() << "some output info" << std::endl;

#endif
}

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


