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

    int N = pow(domainSize / boxSize, DIM);
    int n0 = N / numProc();
    int r = N % numProc();
    if (procID() < r) { n0++; }
    int n = 0;
    for (auto index : layout)
    {
        //pout() << "patch: " << layout[index] << " | local: " << index.local() << " | global: " << index.global() << std::endl;
        n++;
    }
    //pout() << "Iteration Test | N: " << N << " | n: " << n << std::endl;
    EXPECT_EQ(n, n0);
}

#ifdef PR_MPI
TEST(DisjointBoxLayout, LoadBalance) {
    int domainSize = 64;
    int boxSize = 16;
    Box domainBox = Box::Cube(domainSize);
    Point boxSizeVect = Point::Ones(boxSize);
    std::array<bool, DIM> periodicity;
    periodicity.fill(true);
    ProblemDomain domain(domainBox, periodicity);
    DisjointBoxLayout layout_0(domain, boxSizeVect);
    DisjointBoxLayout layout_1(domain, boxSizeVect);
    
    int N = numProc();
    int r = procID();
    int N0 = (N % 2 == 0) ? N / 2 : N / 2 + 1;
    int N1 = N / 2;
    int midRank = N0;
    //std::cout << "N0: " << N0 << " | N1: " << N1 << " | midRank: " << midRank << std::endl;
    layout_0.loadBalance(0, midRank);
    layout_1.loadBalance(midRank, N);

    /*
    for (auto iter : layout_0)
    {
        pout() << "layout_0 | local: " << iter.local() << " | global: " << iter.global() << " | box: " << layout_0[iter] << std::endl;
    }
    for (auto iter : layout_1)
    {
        pout() << "layout_1 | local: " << iter.local() << " | global: " << iter.global() << " | box: " << layout_1[iter] << std::endl;
    }
    */

}
#endif
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
