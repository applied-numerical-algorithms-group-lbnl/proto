#include <gtest/gtest.h>
#include "Proto.H"
#include "Lambdas.H"
#include "BoxOp_MBLaplace.H"
#include "MBMap_XPointRigid.H"

using namespace Proto;

TEST(MBAMR, Construction) {
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

    auto domain = buildXPoint(domainSize);
    std::vector<Point> boxSizeVect(numBlocks, Point::Ones(boxSize));
    std::vector<Point> refRatios(numLevels-1, Point::Ones(refRatio));

    auto coarsePatches = domain.patches(Point::Ones(boxSize));
    MBAMRGrid grid(domain, coarsePatches, boxSizeVect, refRatios);
    for (int li = 1; li < numLevels; li++)
    {
        auto domain = grid[li].domain();
        std::vector<MBPatchID_t> patches;
        Box patchDomain = Box::Cube(domainSize).refine(pow(refRatio,li)).coarsen(boxSize);
        Point patch = patchDomain.high();
        for (int bi = 0; bi < domain.size(); bi++)
        {
            patches.push_back(MBPatchID_t(patch, bi));
        }
        grid[li].define(domain, patches, boxSizeVect);
    }

    MBAMR<BoxOp_MBLaplace, MBMap_XPointRigid, double> amr(
            grid, Point::Ones(refRatio)); 

    MBAMRMap<MBMap_XPointRigid, HOST> map(grid, ghost);
    MBAMRData<double, 1, HOST> data(grid, ghost);
    data.setRandom(0,1);


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
