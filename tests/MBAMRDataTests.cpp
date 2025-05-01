#include <gtest/gtest.h>
#include "Proto.H"
#include "TestFunctions.H"
#include "MBMap_XPointRigid.H"

using namespace Proto;

#if 1
TEST(MBAMRData, Construction) {
    HDF5Handler h5;
    int domainSize = 16;
    int boxSize = 16;
    int numBlocks = 5;
    int numLevels = 3;
    int refRatio = 2;
    int numGhost = 1;
    Array<Point,DIM+1> ghost;
    ghost.fill(Point::Ones(numGhost));

    auto domain = buildXPoint(domainSize);
    std::vector<Point> boxSizeVect(numBlocks, Point::Ones(boxSize));
    std::vector<Point> refRatios(numLevels-1, Point::Ones(refRatio));

    auto coarsePatches = domain.patches(Point::Ones(boxSize));
    MBAMRLayout grid(domain, coarsePatches, boxSizeVect, refRatios);
    for (int li = 1; li < numLevels; li++)
    {
        auto domain = grid[li].domain();
        std::vector<MBPoint> patches;
        Box patchDomain = Box::Cube(domainSize).refine(pow(refRatio,li)).coarsen(boxSize);
        Point patch = patchDomain.high();
        for (int bi = 0; bi < domain.size(); bi++)
        {
            patches.push_back(MBPoint(patch, bi));
        }
        grid[li].define(domain, patches, boxSizeVect);
    }

    MBAMRMap<MBMap_XPointRigid, HOST> map(grid, ghost);
    for (int li = 0; li < numLevels; li++)
    {
        for (BlockIndex bi = 0; bi < numBlocks; bi++)
        {
            map[li][bi].setNumBlocks(numBlocks);
        }
    }
    
    MBAMRData<double, 1, HOST> data(grid, ghost);
    data.setRandom(0,1);
#if PR_VERBOSE > 0
    h5.writeMBAMRData({"data"}, map, data, "MBAMRData_Construction");
#endif

    int numBoxes = domainSize / boxSize * numBlocks;
    for (int li = 0; li < numLevels; li++)
    {
        auto& level = data.getLevel(li);
        auto& layout = grid.getLevel(li);
        for (auto iter : layout)
        {
            auto& patch = level[iter];
            EXPECT_EQ(patch.box(), layout[iter].grow(numGhost));
        }
    }

    for (int li = 0; li < numLevels; li++)
    {
        auto& levelData = data.getLevel(li);
        for (int bi = 0; bi < numBlocks; bi++)
        {
            auto& blockData = data.getBlock(bi);
            auto levelBlockAddress = &(levelData.getBlock(bi));
            auto blockLevelAddress = &(blockData[li]);
            EXPECT_EQ(levelBlockAddress, blockLevelAddress);
        }
    }
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
