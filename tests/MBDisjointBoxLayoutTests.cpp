#include <gtest/gtest.h>
#include "ProtoMMB.H"

#define XPOINT_SIZE 5

using namespace Proto;
MBProblemDomain buildXPoint(int a_domainSize)
{
    MBProblemDomain domain(XPOINT_SIZE);
    auto CCW = CoordPermutation::ccw();
    for (int ii = 0; ii < XPOINT_SIZE; ii++)
    {
        domain.defineBoundary(ii, (ii+1) % XPOINT_SIZE, 0, Side::Hi, CCW);
    }
    for (int bi = 0; bi < XPOINT_SIZE; bi++)
    {
        domain.defineDomain(bi, Point::Ones(a_domainSize));
    }
    return domain;
}

TEST(MBDisjointBoxLayout, Iteration) {
    int domainSize = 64;
    int boxSize = 16;
    auto domain = buildXPoint(domainSize);
    Point boxSizeVect = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);
    
    std::set<std::pair<unsigned int, Box>> correctPatchSet;
    Box patches = Box::Cube(domainSize).coarsen(boxSizeVect);
    for (unsigned int bi = 0; bi < domain.size(); bi++)
    {
        for (auto p : patches)
        {
            Box b(p,p);
            b = b.refine(boxSizeVect);
            std::pair<unsigned int, Box> patch(bi, b);
            correctPatchSet.insert(patch);
        }
    }
    
    int N = pow(domainSize / boxSize, DIM)*domain.size();
    int n0 = N / numProc();
    int r =  N % numProc();
    if (procID() < r) { n0++; }
    std::set<std::pair<unsigned int, Box>> patchSet;
    for (auto index : layout)
    {
        Box b = layout[index];
        auto block = layout.block(index);
        std::pair<unsigned int, Box> patch(block, b);
        EXPECT_EQ(b.sizes(), boxSizeVect);
        EXPECT_EQ(patchSet.find(patch), patchSet.end());
        EXPECT_FALSE(correctPatchSet.find(patch) == correctPatchSet.end());
        for (auto pi : patchSet)
        {
            if (pi.first == block)
            {
                EXPECT_TRUE((pi.second & b).empty());
            }
        }
        patchSet.insert(patch);
    }
    EXPECT_EQ(patchSet.size(), n0);
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
