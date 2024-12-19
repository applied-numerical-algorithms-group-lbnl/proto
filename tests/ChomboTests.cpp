#include <gtest/gtest.h>
#include "Proto.H"
#include "TestFunctions.H"
#define NCOMP 3;

using namespace Proto;

TEST(ChomboTests, GetLayout) {
    int domainSize = 64;
    int boxSize = 16;
    Point boxSizeVect = Point::Ones(boxSize);
    Box domainBox = Box::Cube(domainSize);
    std::array<bool, DIM> periodicity;
    periodicity.fill(true);
    periodicity[0] = false;

    Proto_Transition::DataUtilities<3> utils;
    Box patchDomain = domainBox.coarsen(boxSizeVect);
    std::vector<std::pair<Point, unsigned int>> patches;
    std::map<Point, int> procMap;
    std::vector<int> localSizes(numProc(), 0);
    for (auto pi : patchDomain)
    {
        if (pi == Point::Zeros()) { continue; } //skip this one
        int proc = pi.sum() % numProc();
        patches.push_back(std::make_pair(pi, proc));
        procMap[pi] = proc;
        localSizes[proc]++;
    }
    
    // Build DBL
    ProblemDomain domain(domainBox, periodicity);
    DisjointBoxLayout layout = utils.getLayout(
            patches, boxSizeVect, domain);
    auto sortedLayout = layout.sorted();

    // Test DBL
    EXPECT_EQ(layout.size(), patchDomain.size()-1);
    EXPECT_EQ(sortedLayout.size(), patchDomain.size()-1);
    EXPECT_EQ(layout.localSize(), localSizes[procID()]);
    EXPECT_EQ(sortedLayout.localSize(), localSizes[procID()]);
    for (auto iter : layout)
    {
        Point p = layout.point(iter);
        EXPECT_EQ(procMap[p], procID());
        EXPECT_NE(p, Point::Zeros());
        
    }
    for (auto iter : sortedLayout)
    {
        Point q = sortedLayout.point(iter);
        EXPECT_EQ(procMap[q], procID());
        EXPECT_NE(q, Point::Zeros());
    }
    
    // Exchange Test
    HDF5Handler h5;
    Point ghost = Point::Ones(2);
    LevelBoxData<double, DIM, HOST> data(layout, ghost);
    data.setVal(0);
    for (auto iter : layout)
    {
        auto& patch = data[iter];
        BoxData<double, DIM> tmp(layout[iter]);
        forallInPlace_p(f_pointID, tmp);
        tmp.copyTo(patch);
    }
    data.exchange();
    EXPECT_TRUE(testExchange(data));    

    // Copy Test
    LevelBoxData<double, DIM, HOST> src(layout, ghost);
    LevelBoxData<double, DIM, HOST> dst(layout, ghost);
    LevelBoxData<double, DIM, HOST> sortedSrc(sortedLayout, ghost);
    LevelBoxData<double, DIM, HOST> sortedDst(sortedLayout, ghost);
    src.setVal(0);
    dst.setVal(0);
    for (auto iter : layout)
    {
        auto& patch = src[iter];
        forallInPlace_p(f_pointID, patch);
    }
    for (auto iter : sortedLayout)
    {
        auto& srcPatch = sortedSrc[iter];
        forallInPlace_p(f_pointID, srcPatch);
    }
    src.copyTo(dst);
    src.copyTo(sortedDst);
    
    EXPECT_TRUE(compareLevelData(src, dst));
    EXPECT_TRUE(compareLevelData(sortedSrc, sortedDst));
#if PR_VERBOSE > 0
    h5.writeLevel(1,src,"UNSORTED_SRC");
    h5.writeLevel(1,dst,"UNSORTED_DST");
    h5.writeLevel(1,sortedSrc,"SORTED_SRC");
    h5.writeLevel(1,sortedDst,"SORTED_DST");
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
