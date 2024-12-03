#include <gtest/gtest.h>
#include "Proto.H"
#include "TestFunctions.H"
#include "MBMap_XPointRigid.H"

using namespace Proto;

namespace {
    
    std::map<Point, Box> telescopingCFBoundary_Coarse(int domainSize)
    {
        std::map<Point, Box> boundMap;
        Box domainBox = Box::Cube(domainSize);

        Point refinedRegionLow = Point::Ones(domainSize / 2);
        Point refinedRegionHigh = Point::Ones(domainSize - 1);
        Box refinedRegion(refinedRegionLow, refinedRegionHigh);

        for (auto dir : Point::DirectionsOfCodim(1))
        {
            boundMap[dir] = Box();
        }
        boundMap[Point::X()] = refinedRegion.adjacent(0, Side::Lo, 1);
        boundMap[Point::Y()] = refinedRegion.adjacent(1, Side::Lo, 1);

        return boundMap;
    }

    std::map<Point, Box> telescopingCFBoundary_Fine(int domainSize)
    {
        std::map<Point, Box> boundMap;
        Box domainBox = Box::Cube(domainSize);

        Point refinedRegionLow = Point::Ones(domainSize / 2);
        Point refinedRegionHigh = Point::Ones(domainSize - 1);
        Box refinedRegion(refinedRegionLow, refinedRegionHigh);

        for (auto dir : Point::DirectionsOfCodim(1))
        {
            boundMap[dir] = Box();
        }
        boundMap[-Point::X()] = refinedRegion.adjacent(0, Side::Lo, 1);
        boundMap[-Point::Y()] = refinedRegion.adjacent(1, Side::Lo, 1);

        return boundMap;
    }
    
    std::map<Point, Box> refinedBlockBoundaryCFBoundary_Coarse(int domainSize, int boxSize, int refRatio)
    {
        std::map<Point, Box> boundMap;
        Box domainBox = Box::Cube(domainSize);

        int refinedRegionWidth = boxSize / refRatio;
        Box unrefinedRegion = Box::Cube(domainSize).extrude(Point::X(), -refinedRegionWidth);

        for (auto dir : Point::DirectionsOfCodim(1))
        {
            boundMap[dir] = Box();
        }
        boundMap[Point::X()] = unrefinedRegion.edge(0, Side::Hi, 1);
        boundMap[Point::Y()] = unrefinedRegion.edge(1, Side::Hi, 1);

        return boundMap;
    }
    
    std::map<Point, Box> refinedBlockBoundaryCFBoundary_Fine(int domainSize, int boxSize, int refRatio)
    {
        std::map<Point, Box> boundMap;
        Box domainBox = Box::Cube(domainSize);

        int refinedRegionWidth = boxSize / refRatio;
        Box unrefinedRegion = Box::Cube(domainSize).extrude(Point::X(), -refinedRegionWidth);
        Box refinedRegion = unrefinedRegion.adjacent(0, Side::Hi, refinedRegionWidth);

        for (auto dir : Point::DirectionsOfCodim(1))
        {
            boundMap[dir] = Box();
        }
        boundMap[Point::X()]  = refinedRegion.adjacent(0, Side::Hi, 1).extrude(Point::Y(), -refinedRegionWidth);
        boundMap[-Point::X()] = unrefinedRegion.edge(0, Side::Hi, 1);

        return boundMap;
    }

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
        counterData.printData();
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
        h5.writeMBAMRData({"data"}, map, data, "Telescoping_XPoint_Data");
    #endif

    Array<double, DIM> gridSpacing = Point::Ones();
    gridSpacing /= domainSize;
    MBLevelFluxRegister<double, 1, HOST> fluxRegister(grid[0], grid[1], refRatios, gridSpacing);
    MBLevelFluxRegisterTester<double, 1, HOST> tester(fluxRegister);
    
    auto cfBounds_Coarse = telescopingCFBoundary_Coarse(domainSize);

    for (auto iter : grid[0])
    {
        auto coarseRegisters = tester.getCoarseRegistersAtIndex(iter);
        for (auto dir : Point::DirectionsOfCodim(1))
        {
        EXPECT_TRUE(checkRegisters(coarseRegisters, cfBounds_Coarse[dir], dir));
        }
    }

    auto cfBounds_Fine = telescopingCFBoundary_Fine(domainSize);
    for (auto iter : grid[1])
    {
        PatchID patch = grid[1].point(iter);
        auto fineRegisters = tester.getFineRegistersAtIndex(iter);

        for (auto dir : Point::DirectionsOfCodim(1))
        {
            Box cfBound = grid[1][iter].adjacent(dir).coarsen(refRatio);
            cfBound &= cfBounds_Fine[dir];
            EXPECT_TRUE(checkRegisters(fineRegisters, cfBound, dir));
        }
    }
}

TEST(MBLevelFluxRegister, RefinedBlockBoundaryXPointConstruction) {
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

    auto grid = refinedBlockBoundaryXPointGrid(domainSize, numLevels, refRatio, boxSize);
    MBAMRMap<MBMap_XPointRigid, HOST> map(grid, ghost);
    MBAMRData<double, 1, HOST> data(grid, ghost);

    #if PR_VERBOSE > 0
        h5.writeMBAMRData({"data"}, map, data, "BlockBoundary_XPoint_Data");
    #endif

    Array<double, DIM> gridSpacing = Point::Ones();
    gridSpacing /= domainSize;
    MBLevelFluxRegister<double, 1, HOST> fluxRegister(grid[0], grid[1], refRatios, gridSpacing);
    MBLevelFluxRegisterTester<double, 1, HOST> tester(fluxRegister);
    
    auto cfBounds_Coarse = refinedBlockBoundaryCFBoundary_Coarse(domainSize, boxSize, refRatio);
    for (auto iter : grid[0])
    {
        auto coarseRegisters = tester.getCoarseRegistersAtIndex(iter);
        for (auto dir : Point::DirectionsOfCodim(1))
        {
            EXPECT_TRUE(checkRegisters(coarseRegisters, cfBounds_Coarse[dir], dir));
        }
    }

    auto cfBounds_Fine = refinedBlockBoundaryCFBoundary_Fine(domainSize, boxSize, refRatio);
    for (auto iter : grid[1])
    {
        PatchID patch = grid[1].point(iter);
        auto fineRegisters = tester.getFineRegistersAtIndex(iter);

        for (auto dir : Point::DirectionsOfCodim(1))
        {
            Box cfBound = grid[1][iter].adjacent(dir).coarsen(refRatio);
            cfBound &= cfBounds_Fine[dir];
            EXPECT_TRUE(checkRegisters(fineRegisters, cfBound, dir));
        }
    }
}

TEST(MBLevelFluxRegister, TelescopingXPointIncrement) {
    int domainSize = 16;
    int boxSize = 16;
    int refRatio = 4;
    int ghostWidth = 2;

    Point refRatios = Point::Ones(refRatio);
    Point ghostWidths = Point::Ones(ghostWidth);

    auto grid = telescopingXPointGrid(domainSize, 2, refRatio, boxSize);
    MBAMRMap<MBMap_XPointRigid, HOST> map(grid, ghostWidths);
    MBAMRData<double, 1, HOST> data(grid, ghostWidths);
}

#if DIM == 3
TEST(MBLevelFluxRegister, Reflux_CubedSphere) {
   
    int domainSize = 32;
    int boxSize = 16;
    int ghostSize = 2;

}
#endif
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
