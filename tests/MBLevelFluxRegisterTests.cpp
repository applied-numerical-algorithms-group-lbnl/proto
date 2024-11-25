#include <gtest/gtest.h>
#include "Proto.H"
#include "TestFunctions.H"
#include "MBMap_XPointRigid.H"

using namespace Proto;

namespace {
    template<typename T, unsigned int C, MemType MEM>
    bool checkRegisters(
        std::vector<Register<T,C,MEM>>& a_registers,
        Box a_coarseFineBoundary,
        Point a_dir)
    {
        BoxData<T,C,MEM> counterData(a_coarseFineBoundary);
        counterData.setVal(0);
        
        for (auto ri : a_registers)
        {
            if (ri.dir() == a_dir)
            {
                BoxData<T,C,MEM> rdata(ri.data().box());
                if (!a_coarseFineBoundary.containsBox(rdata.box())) { return false; }
                rdata.setVal(1);
                counterData += rdata;
            }
        }
        int numDataInCoarseFineBoundary = a_coarseFineBoundary.size()*C;
        bool success = (counterData.sum() == numDataInCoarseFineBoundary);
        return success;
    }

}

#if PR_MMB
TEST(MBLevelFluxRegister, TelescopingXPointConstruction) {
    #if PR_VERBOSE > 0
        HDF5Handler h5;
    #endif

    int domainSize = 16;
    int boxSize = domainSize;
    int numBlocks = MB_MAP_XPOINT_NUM_BLOCKS;
    int numLevels = 2;
    int refRatio = 4;
    int numGhost = 2;
    
    Point refRatios = Point::Ones(refRatio);

    Array<Point,DIM+1> ghost;
    ghost.fill(Point::Ones(numGhost));

    auto grid = telescopingXPointGrid(domainSize, numLevels, refRatio, boxSize);
    MBAMRMap<MBMap_XPointRigid, HOST> map(grid, ghost);
    MBAMRData<double, 1, HOST> data(grid, ghost);

    #if PR_VERBOSE > 0
        h5.writeMBAMRData({"data"}, map, data, "Reflux_XPoint_Data");
    #endif

    Array<double, DIM> gridSpacing = Point::Ones();
    gridSpacing /= domainSize;
    MBLevelFluxRegister<double, 1, HOST> fluxRegister(grid[0], grid[1], refRatios, gridSpacing);
    MBLevelFluxRegisterTester<double, 1, HOST> tester(fluxRegister);

    int fineDomainSize = domainSize * refRatio;
    int finePatchesPerBoundary = (fineDomainSize / 2) / boxSize; 
    Point refinedRegionLow = Point::Ones(domainSize / 2);
    Point refinedRegionHigh = Point::Ones(domainSize - 1);
    Box refinedRegion(refinedRegionLow, refinedRegionHigh);
    Box coarseFineBound_X = refinedRegion.adjacent(0, Side::Lo, 1);
    Box coarseFineBound_Y = refinedRegion.adjacent(1, Side::Lo, 1);
    Box emptyBox;
    

    for (auto iter : grid[0])
    {
        auto coarseRegisters = tester.getCoarseRegistersAtIndex(iter);
        EXPECT_EQ(coarseRegisters.size(), 2*finePatchesPerBoundary);
        EXPECT_TRUE(checkRegisters(coarseRegisters, coarseFineBound_X, Point::X()));
        EXPECT_TRUE(checkRegisters(coarseRegisters, coarseFineBound_Y, Point::Y()));
        EXPECT_TRUE(checkRegisters(coarseRegisters, emptyBox, -Point::X()));
        EXPECT_TRUE(checkRegisters(coarseRegisters, emptyBox, -Point::Y()));
    }

    Box refinedPatches = refinedRegion.coarsen(boxSize);
    for (auto iter : grid[1])
    {
        PatchID patch = grid[1].point(iter);

        Box cfBound_x = grid[1][iter].adjacent(0, Side::Lo, 1).coarsen(refRatio);
        Box cfBound_y = grid[1][iter].adjacent(1, Side::Lo, 1).coarsen(refRatio);
        cfBound_x &= coarseFineBound_X;
        cfBound_y &= coarseFineBound_Y;
        auto fineRegisters = tester.getFineRegistersAtIndex(iter);
        //EXPECT_EQ(fineRegisters.size(), 2*finePatchesPerBoundary);
        EXPECT_TRUE(checkRegisters(fineRegisters, cfBound_x, -Point::X()));
        EXPECT_TRUE(checkRegisters(fineRegisters, cfBound_y, -Point::Y()));
        EXPECT_TRUE(checkRegisters(fineRegisters, emptyBox, Point::X()));
        EXPECT_TRUE(checkRegisters(fineRegisters, emptyBox, Point::Y()));
    }
    
}
#endif

#if DIM == 3
TEST(MBLevelFluxRegister, Reflux_CubedSphere) {
   
    int domainSize = 32;
    int boxSize = 16;
    int ghostSize = 2;

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
