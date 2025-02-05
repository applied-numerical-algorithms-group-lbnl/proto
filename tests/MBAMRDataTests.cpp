#include <gtest/gtest.h>
#include "Proto.H"
#include "TestFunctions.H"
#include "MBMap_XPointRigid.H"

using namespace Proto;

namespace {
    MBAMRGrid testGrid(int domainSize, int boxSize, int refRatio, int numLevels, int numBlocks)
    {
        auto domain = buildXPoint(domainSize, numBlocks);
        std::vector<Point> boxSizeVect(domain.numBlocks(), Point::Ones(boxSize));
        std::vector<Point> refRatios(numLevels-1, Point::Ones(refRatio));

        auto coarsePatches = domain.patches(Point::Ones(boxSize));
        MBAMRGrid grid(domain, coarsePatches, boxSizeVect, refRatios);
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
        return grid;
    }
}

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
#if 1
TEST(MBAMRData, Interpolate) {

    HDF5Handler h5;
    int domainSize = 16;
    int boxSize = 16;
    int numBlocks = 5;
    int numLevels = 3;
    int refRatio = 2;
    int numGhost = 4;
    //Array<Point,DIM+1> ghost;
    //ghost.fill(Point::Ones(numGhost));
    Point ghost = Point::Ones(numGhost);
    auto grid = testGrid(domainSize, boxSize, refRatio, numLevels, numBlocks);
    MBAMRMap<MBMap_XPointRigid, HOST> map(grid, ghost);
    pr_out() << "num ghost: " << ghost << std::endl;
    MBAMRData<double, 1, HOST> data(grid, ghost);
    Array<double, DIM> X0;
    X0.fill(0);
    
    for (auto iter : grid[0])
    {
        auto& patch = data.getLevel(0)[iter];
        auto& X = map.getLevel(0).map()[iter];
        forallInPlace(f_bell, patch, X, X0);
    }
    for (int lvl = 1; lvl < numLevels; lvl++)
    {
        data.getLevel(lvl).setVal(-1);
    }
    
    h5.writeMBAMRData({"phi"}, map, data, "MBAMRDataTests_Interpolate_DATA_0");
    data.interpolate();
    h5.writeMBAMRData({"phi"}, map, data, "MBAMRDataTests_Interpolate_DATA_1");
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
