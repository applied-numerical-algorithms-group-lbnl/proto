#include <gtest/gtest.h>
#include "Proto.H"
#include "Lambdas.H"

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
    ghost[0] = Point::Ones(numGhost);

    auto domain = buildXPoint(domainSize, numBlocks);
    std::vector<Point> boxSizeVect(numBlocks, Point::Ones(boxSize));
    std::vector<Point> refRatios(numLevels-1, Point::Ones(refRatio));

    auto coarsePatches = domain.patches(Point::Ones(boxSize));
    MBAMRGrid grid(domain, coarsePatches, boxSizeVect, refRatios);
    MBAMRData<double, 1, HOST> data(grid, ghost);
    data.setRandom(0,1);
#if PR_VERBOSE > 0
    h5.writeMBAMRData({"data"}, data, "DATA");
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
