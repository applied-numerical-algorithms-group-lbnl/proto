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
        int numDataInCoarseFineBoundary = a_coarseFineBoundary.size()*C;
        bool success = (counterData.sum() == numDataInCoarseFineBoundary);
        return success;
    }


    template<typename T, unsigned int C, MemType MEM>
    PROTO_KERNEL_START
    void f_testFluxF(Point& a_pt, Var<T,C,MEM>& a_F, int a_dir)
    {
        double sign = 1.0;// pow(-1,a_dir);
        for (int cc = 0; cc < C; cc++)
        {
            a_F(cc) = sign * (a_pt.sum()*1.0 + 1000.0*cc);
        }

    }
    PROTO_KERNEL_END(f_testFluxF, f_testFlux);

    template<typename T, unsigned int C>
    BoxData<T,C> testCoarseReflux(Box registerBox, Point dir)
    {
        PROTO_ASSERT(dir.codim() == 1, "Error");
        int coord = dir.firstNonZeroIndex();
        auto side = (dir[coord] > 0) ? Side::Hi : Side::Lo;
        Stencil<T> shiftAndScale;
        if (side == Side::Hi)
        {
            // Hi side flux is outward flowing -> contributes negatively to divergence
            shiftAndScale = -1.0*Shift::Basis(coord,+1);
        } else {
            // Lo side flux is inward flowing -> contributes positively to divergence
            shiftAndScale = +1.0*Shift::Zeros();
        }

        // Coarse flux contribution is subtracted during refluxing
        shiftAndScale *= -1;

        Box sourceBox = registerBox.grow(PR_NODE);
        auto sourceData = forall_p<T,C>(f_testFlux, sourceBox, coord);
        BoxData<T,C> outData(registerBox);
        Box rangeBox = shiftAndScale.range(sourceBox);
        PROTO_ASSERT(rangeBox.containsBox(registerBox), "Error");
        outData |= shiftAndScale(sourceData);
        return outData;
    }

    template<typename T, unsigned int C>
    BoxData<T,C> testFineReflux(Box registerBox, Point dir, int refRatio)
    {
        PROTO_ASSERT(dir.codim() == 1, "Error");
        int coord = dir.firstNonZeroIndex();
        auto side = (dir[coord] > 0) ? Side::Hi : Side::Lo;
        auto avgShiftAndScale = Stencil<T>::AvgDownFace(coord, Side::Lo, refRatio);
        if (side == Side::Hi)
        {
            
        } else {
            // low side fluxes contribute negatively to the divergence of the adjacent cell
            avgShiftAndScale.destShift() = Point::Basis(coord,Side::Lo);
            avgShiftAndScale *= -1;
        }

        Box sourceBox = registerBox.grow(PR_NODE);
        sourceBox = sourceBox.refine(refRatio);
        auto sourceData = forall_p<T,C>(f_testFlux, sourceBox, coord);
        BoxData<T,C> outData(registerBox);
        Box rangeBox = avgShiftAndScale.range(sourceBox);
        PROTO_ASSERT(rangeBox.containsBox(registerBox), "Error");
        outData |= avgShiftAndScale(sourceData);
        return outData;
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
    constexpr int NUM_COMPS = 2;
    
    int domainSize = 16;
    int boxSize = 16;
    int refRatio = 4;
    int ghostWidth = 2;

    Point refRatios = Point::Ones(refRatio);
    Point ghostWidths = Point::Ones(ghostWidth);

    auto grid = telescopingXPointGrid(domainSize, 2, refRatio, boxSize);

    Array<double, DIM> gridSpacing = Point::Ones();
    gridSpacing /= domainSize;
    MBLevelFluxRegister<double, NUM_COMPS, HOST> coarseFluxRegister(grid[0], grid[1], refRatios, gridSpacing);
    MBLevelFluxRegister<double, NUM_COMPS, HOST> fineFluxRegister(grid[0], grid[1], refRatios, gridSpacing);
    MBLevelFluxRegister<double, NUM_COMPS, HOST> refluxFluxRegister(grid[0], grid[1], refRatios, gridSpacing);

    for (auto iter : grid[0])
    {
        FluxBoxData<double, NUM_COMPS> fluxes(grid[0][iter]);
        for (int dd = 0; dd < DIM; dd++)
        {
            forallInPlace_p(f_testFlux, fluxes[dd], dd);
            coarseFluxRegister.incrementCoarseRegister(fluxes[dd], iter, dd);
            refluxFluxRegister.incrementCoarseRegister(fluxes[dd], iter, dd);
        }
    }

    for (auto iter : grid[1])
    {
        FluxBoxData<double, NUM_COMPS> fluxes(grid[1][iter]);
        for (int dd = 0; dd < DIM; dd++)
        {
            forallInPlace_p(f_testFlux, fluxes[dd], dd);
            fineFluxRegister.incrementFineRegister(fluxes[dd], iter, dd);
            refluxFluxRegister.incrementFineRegister(fluxes[dd], iter, dd);
        }
    }

    MBLevelBoxData<double, NUM_COMPS, HOST> refluxRegisters(grid[0], Point::Zeros());
    MBLevelBoxData<double, NUM_COMPS, HOST> coarseRegisters(grid[0], Point::Zeros());
    MBLevelBoxData<double, NUM_COMPS, HOST> fineRegisters(grid[0], Point::Zeros());

    coarseRegisters.setVal(0);
    fineRegisters.setVal(0);
    refluxRegisters.setVal(0);

    coarseFluxRegister.applyRefluxCorrection(coarseRegisters, 1.0);
    fineFluxRegister.applyRefluxCorrection(fineRegisters, 1.0);
    refluxFluxRegister.applyRefluxCorrection(refluxRegisters, 1.0);

#if PR_VERBOSE > 0
    HDF5Handler h5;

    MBAMRMap<MBMap_XPointRigid, HOST> map(grid, ghostWidths);

    h5.writeMBLevel(map[0], coarseRegisters, "COARSE_REGISTERS");
    h5.writeMBLevel(map[0], fineRegisters, "FINE_REGISTERS");
    h5.writeMBLevel(map[0], refluxRegisters, "REFLUX_REGISTERS");
#endif
    auto coarseCFBounds = telescopingCFBoundary_Coarse(domainSize);
    for (auto iter : grid[0])
    {
        Box B0 = grid[0][iter];
        auto& coarseData = coarseRegisters[iter];
        BoxData<double, NUM_COMPS> error(coarseData.box());
        error.setVal(0);
        error += coarseData; 
        for (auto dir : Point::DirectionsOfCodim(1))
        {
            Box registerBox = coarseCFBounds[dir];
            if (registerBox.empty()) { continue; }
            auto coarseSoln = testCoarseReflux<double, NUM_COMPS>(registerBox, dir);
            error -= coarseSoln; 
        }
        
        EXPECT_LT(error.absMax(), 1e-12);
    }

    auto fineCFBounds = telescopingCFBoundary_Fine(domainSize);
    for (auto iter : grid[0])
    {
        Box B0 = grid[0][iter];
        auto& fineData = fineRegisters[iter];
        BoxData<double, NUM_COMPS> error(fineData.box());
        error.setVal(0);
        error += fineData; 
        for (auto dir : Point::DirectionsOfCodim(1))
        {
            Box registerBox = fineCFBounds[dir];
            if (registerBox.empty()) { continue; }
            auto fineSoln = testFineReflux<double, NUM_COMPS>(registerBox, dir, refRatio);
            error -= fineSoln; 
        }
        
        EXPECT_LT(error.absMax(), 1e-12);
    }

    for (auto iter : grid[0])
    {
        Box B0 = grid[0][iter];
        auto& refluxData = refluxRegisters[iter];
        BoxData<double, NUM_COMPS> error(refluxData.box());
        error.setVal(0);
        error += refluxData; 
        for (auto dir : Point::DirectionsOfCodim(1))
        {
            Box registerBox = coarseCFBounds[dir];
            if (registerBox.empty()) { continue; }
            BoxData<double, NUM_COMPS> refluxSoln(registerBox);
            refluxSoln.setVal(0);
            auto coarseSoln = testCoarseReflux<double, NUM_COMPS>(registerBox, dir);
            auto fineSoln = testFineReflux<double, NUM_COMPS>(registerBox, -dir, refRatio);
            
            refluxSoln += fineSoln; 
            refluxSoln += coarseSoln;

            error -= refluxSoln;
        }
        
        EXPECT_LT(error.absMax(), 1e-12);
    }

}

#if DIM == 3
TEST(MBLevelFluxRegister, Reflux_CubedSphere) {
   
    int domainSize = 32;
    int boxSize = 16;
    int ghostWidth = 2;

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
