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

    std::set<Box> correctBoxSet;
    Box patches = domainBox.coarsen(boxSizeVect);
    for (auto p : patches)
    {
        Box b(p,p);
        correctBoxSet.insert(b.refine(boxSizeVect));
    }

    int N = pow(domainSize / boxSize, DIM);
    int n0 = N / numProc();
    int r = N % numProc();
    if (procID() < r) { n0++; }
    int n = 0;
    std::set<Box> boxSet;
    for (auto index : layout)
    {
        n++;
        Box b = layout[index];
        EXPECT_EQ(b.sizes(), boxSizeVect);
        EXPECT_EQ(boxSet.find(b), boxSet.end());
        EXPECT_FALSE(correctBoxSet.find(b) == correctBoxSet.end());
        for (auto bi : boxSet)
        {
            EXPECT_TRUE((bi & b).empty());
        }
        boxSet.insert(b);
    }
    EXPECT_EQ(n, n0);
}

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
    
    std::set<Box> correctBoxSet;
    Box patches = domainBox.coarsen(boxSizeVect);
    for (auto p : patches)
    {
        Box b(p,p);
        correctBoxSet.insert(b.refine(boxSizeVect));
    }

    int N = numProc();
    int r = procID();
    int N0 = (N % 2 == 0) ? N / 2 : N / 2 + 1;
    int N1 = N / 2;
    int midRank = N0;
    layout_0.loadBalance(0, midRank);
    layout_1.loadBalance(midRank, N);
   
    int numBoxes = pow(domainSize / boxSize, DIM);
    int num_0 = (numBoxes / midRank);
    int rem_0 = (numBoxes % midRank);
    num_0 = procID() < rem_0 ? num_0+1 : num_0;
    if (procID() >= midRank) { num_0 = 0; }

    int num_1 = 0;
    if (N > 1)
    {
        num_1 = (numBoxes / (N - midRank));
        int rem_1 = (numBoxes % (N - midRank));
        num_1 = procID() < rem_1 ? num_1+1 : num_1;
    }
    if (procID() < midRank) { num_1 = 0; }

    int n0 = 0;
    int n1 = 0;
    std::set<Box> boxSet_0;
    std::set<Box> boxSet_1;
    for (auto iter : layout_0)
    {
        n0++;
        Box b = layout_0[iter];
        EXPECT_EQ(b.sizes(), boxSizeVect);
        EXPECT_EQ(boxSet_0.find(b), boxSet_0.end());
        EXPECT_FALSE(correctBoxSet.find(b) == correctBoxSet.end());
        for (auto bi : boxSet_0)
        {
            EXPECT_TRUE((bi & b).empty());
        }
        boxSet_0.insert(b);
    }
    for (auto iter : layout_1)
    {
        n1++;
        Box b = layout_1[iter];
        EXPECT_EQ(b.sizes(), boxSizeVect);
        EXPECT_EQ(boxSet_1.find(b), boxSet_1.end());
        EXPECT_FALSE(correctBoxSet.find(b) == correctBoxSet.end());
        for (auto bi : boxSet_0)
        {
            EXPECT_TRUE((bi & b).empty());
        }
        boxSet_1.insert(b);
    }
    EXPECT_EQ(n0, num_0);
    EXPECT_EQ(n1, num_1);
}

TEST(DisjointBoxLayout, LoadAssign) {
    int domainSize = 64;
    int boxSize = 16;
    Box domainBox = Box::Cube(domainSize);
    Point boxSizeVect = Point::Ones(boxSize);
    std::array<bool, DIM> periodicity;
    periodicity.fill(true);
    ProblemDomain domain(domainBox, periodicity);
    DisjointBoxLayout layout(domain, boxSizeVect);
    
    std::set<Box> correctBoxSet;
    Box patches = domainBox.coarsen(boxSizeVect);
    for (auto p : patches)
    {
        Box b(p,p);
        correctBoxSet.insert(b.refine(boxSizeVect));
    }

    unsigned int numBoxes = pow(domainSize / boxSize, DIM);
    unsigned int numUnassigned = numBoxes;
    std::vector<std::pair<int, unsigned int>> procAssign;
    procAssign.resize(numProc());
    for (int ii = 0; ii < numProc(); ii++)
    {
        unsigned int assignNum = min(ii+1, (int)numUnassigned);
        procAssign[ii] = std::pair<int, unsigned int>(ii, assignNum);
        numUnassigned -= assignNum;
    }
    procAssign[numProc()-1].second += numUnassigned;
    layout.loadAssign(procAssign);
    int n0 = 0;
    for (auto data : procAssign)
    {
        if (data.first == procID()) { n0 = data.second; }
    }

    int n = 0;
    std::set<Box> boxSet;
    for (auto ind : layout)
    {
        Box b = layout[ind];
        EXPECT_EQ(b.sizes(), boxSizeVect);
        EXPECT_EQ(boxSet.find(b), boxSet.end());
        EXPECT_FALSE(correctBoxSet.find(b) == correctBoxSet.end());
        for (auto bi : boxSet)
        {
            EXPECT_TRUE((bi & b).empty());
        }
        boxSet.insert(b);
        n++;
    }
    EXPECT_EQ(n, n0);
}

TEST(DisjointBoxLayout, Contains) {
    int domainSize = 64;
    int boxSize = 16;
    Box domainBox = Box::Cube(domainSize);
    Point boxSizeVect = Point::Ones(boxSize);
    std::array<bool, DIM> periodicity;
    periodicity.fill(false);
    periodicity[0] = true; //periodic in X
    ProblemDomain domain(domainBox, periodicity);
    DisjointBoxLayout layout(domain, boxSizeVect);

    Box testPatches = domainBox.coarsen(boxSizeVect).grow(1);
    Box validPatches = domainBox.coarsen(boxSizeVect).grow(0,1);

    for (auto p : testPatches)
    {
        if (validPatches.contains(p))
        {
            EXPECT_TRUE(layout.contains(p));
        } else {
            EXPECT_FALSE(layout.contains(p));
        }
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
