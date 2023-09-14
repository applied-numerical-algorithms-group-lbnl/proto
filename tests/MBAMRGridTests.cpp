#include <gtest/gtest.h>
#include "Proto.H"
#include "Lambdas.H"

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
        auto& layout = grid.getLevel(li);
        EXPECT_EQ(layout.numBoxes(), numBoxes);
        for (auto iter : layout)
        {
            auto block = layout.block(iter);
            auto box = layout[iter];
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
