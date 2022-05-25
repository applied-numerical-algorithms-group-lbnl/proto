#include <gtest/gtest.h>
#include "Proto.H"

using namespace Proto;

TEST(DisjointBoxLayout, Iteration) {
    int domainSize = 64;
    int boxSize = 16;
    Box domainBox = Box::Cube(domainSize);
    Point boxSizeVect = Point::Ones(boxSize);
    std::array<bool, DIM> periodicity;
    periodicity.fill(true);
    ProblemDomain domain(domainBox, periodicity);
    DisjointBoxLayout layout(domain, boxSizeVect);

    for (auto index : layout)
    {
        pout() << layout[index] << std::endl;
    }
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
