#include <gtest/gtest.h>
#include "Proto.H"
#include "TestFunctions.H"

using namespace Proto;

#if 1
TEST(MBAMRGrid, Construction) {
    int domainSize = 16;
    int boxSize = 16;
    int numBlocks = 5;
    int numLevels = 3;
    int refRatio = 2;
    auto domain = buildXPoint(domainSize, numBlocks);
    std::vector<Point> boxSizeVect(numBlocks, Point::Ones(boxSize));
    std::vector<Point> refRatios(numLevels-1, Point::Ones(refRatio));

    auto coarsePatches = domain.patches(Point::Ones(boxSize));
    MBAMRGrid grid(domain, coarsePatches, boxSizeVect, refRatios);
    int numBoxes = domainSize / boxSize * numBlocks;
    for (int li = 0; li < numLevels; li++)
    {
        auto& levelLayout = grid.getLevel(li);
        EXPECT_EQ(levelLayout.numBoxes(), numBoxes);
        for (int bi = 0; bi < numBlocks; bi++)
        {
            auto& blockGrid = grid.getBlock(bi); 
            EXPECT_TRUE(levelLayout.getBlock(bi) == blockGrid[li]);
            EXPECT_EQ(levelLayout.getBlock(bi).partition(), blockGrid[li].partition());
        }

        for (auto iter : levelLayout)
        {
            auto block = levelLayout.block(iter);
            auto box = levelLayout[iter];
            EXPECT_EQ(box.sizes(), boxSizeVect[block]);
        }
        numBoxes *= pow(refRatio, DIM);
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
