#include <gtest/gtest.h>
#include "Proto.H"
#include "Lambdas.H"

using namespace Proto;

TEST(MBMultigridTests, Construction) {
    
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
