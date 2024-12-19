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

        Point refinedRegionHigh = Point::Ones(domainSize - 1);
        Point refinedRegionLow = Point::Zeros();
        refinedRegionLow[0] = domainSize / 2;
        refinedRegionLow[1] = domainSize / 2;
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

        Point refinedRegionHigh = Point::Ones(domainSize - 1);
        Point refinedRegionLow = Point::Zeros();
        refinedRegionLow[0] = domainSize / 2;
        refinedRegionLow[1] = domainSize / 2;
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
                if (!a_coarseFineBoundary.containsBox(rdata.box()))
                {
                    return false;
                }
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
        for (int cc = 0; cc < C; cc++)
        {
            a_F(cc) = a_pt[1] + 100*a_pt[0] + 10000*cc;
        }

    }
    PROTO_KERNEL_END(f_testFluxF, f_testFlux);

    template<typename T, unsigned int C>
    BoxData<T,C> testCoarseReflux(Box registerBox, Point dir, Array<T,DIM> gridSpacing)
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
        shiftAndScale *= (1.0/gridSpacing[coord]);

        Box sourceBox = registerBox.grow(PR_NODE);
        auto sourceData = forall_p<T,C>(f_testFlux, sourceBox, coord);
        BoxData<T,C> outData(registerBox);
        Box rangeBox = shiftAndScale.range(sourceBox);
        PROTO_ASSERT(rangeBox.containsBox(registerBox), "Error");
        outData |= shiftAndScale(sourceData);
        return outData;
    }

    template<typename T, unsigned int C>
    BoxData<T,C> testFineReflux(Box registerBox, Point dir, int refRatio, Array<T,DIM> gridSpacing)
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

        avgShiftAndScale *= (1.0/gridSpacing[coord]);

        Box sourceBox = registerBox.grow(PR_NODE);
        sourceBox = sourceBox.refine(refRatio);
        auto sourceData = forall_p<T,C>(f_testFlux, sourceBox, coord);
        BoxData<T,C> outData(registerBox);
        Box rangeBox = avgShiftAndScale.range(sourceBox);
        PROTO_ASSERT(rangeBox.containsBox(registerBox), "Error");
        outData |= avgShiftAndScale(sourceData);
        return outData;
    }

    template<typename T, unsigned int C, MemType MEM>
    MBIndex rotateFluxFromFineRegister(
        BoxData<T,C,MEM>& flux,
        MBIndex fineIndex,
        Point dir,
        MBDisjointBoxLayout fineLayout)
    {
        PatchID finePatch = fineLayout.point(fineIndex);
        BlockIndex block = fineLayout.block(fineIndex);
        bool patchOnBlockBoundary = fineLayout.onBlockBoundary(finePatch, block, dir);
        if (patchOnBlockBoundary)
        {
            BlockIndex adjBlock = fineLayout.domain().graph().adjacent(block, dir);
            CoordPermutation R = fineLayout.domain().graph().rotation(block, dir);
            Box rotatedDomain = fineLayout.domain().convertBox(flux.box(), block, adjBlock);
            flux.rotate(rotatedDomain, R);

            PatchID adjPatch = rotatedDomain.coarsen(fineLayout.boxSizes()[adjBlock]).low();
            auto adjIndex = fineLayout.find(adjPatch, adjBlock);
            if (adjIndex == *fineLayout.end())
            {
                std::cout << "Couldn't find adjacent index" << std::endl;
                std::cout << "source data: " << finePatch << " | " << block << std::endl;
                std::cout << "dir: " << dir << std::endl;
                std::cout << "adj data: " << adjPatch << " | " << adjBlock << std::endl;
            }
        }
        return fineLayout.find(finePatch, block);
    }

    std::tuple<Point, Box, BlockIndex> getFineRegisterArgsFromCoarse(
        Point coarseDir,
        Box coarseRegisterBox,
        MBDisjointBoxLayout coarseLayout,
        MBIndex coarseIndex)
    {
        auto block = coarseLayout.block(coarseIndex);
        Box domainBox = coarseLayout.domain().getBlock(block).box();
        bool boxOnBlockBoundary = domainBox.edge(coarseDir).containsBox(coarseRegisterBox);
        if (boxOnBlockBoundary)
        {
            auto adjBlock = coarseLayout.domain().graph().adjacent(block, coarseDir);
            Box fineRegisterBox = coarseLayout.domain().convertBox(coarseRegisterBox, block, adjBlock);
            Point fineDir = coarseLayout.domain().graph().reverseDir(block, adjBlock, coarseDir);
            return std::make_tuple(fineDir, fineRegisterBox, adjBlock);
        } else {
            return std::make_tuple(-coarseDir, coarseRegisterBox, block);
        }
    }

    template <typename T, unsigned int C, MemType MEM>
    std::tuple<bool, bool, bool> initializeFluxRegisters(
        MBLevelFluxRegister<T,C,MEM>& coarseFluxRegister,
        MBLevelFluxRegister<T,C,MEM>& fineFluxRegister,
        MBLevelFluxRegister<T,C,MEM>& refluxFluxRegister,
        MBAMRGrid& grid)
    {
        for (auto iter : grid[0])
        {
            FluxBoxData<T, C> fluxes(grid[0][iter]);
            for (int dd = 0; dd < DIM; dd++)
            {
                forallInPlace_p(f_testFlux, fluxes[dd], dd);
                coarseFluxRegister.incrementCoarseRegister(fluxes[dd], iter, dd);
                refluxFluxRegister.incrementCoarseRegister(fluxes[dd], iter, dd);
            }
        }

        for (auto iter : grid[1])
        {
            FluxBoxData<T,C> fluxes(grid[1][iter]);
            for (int dd = 0; dd < DIM; dd++)
            {
                forallInPlace_p(f_testFlux, fluxes[dd], dd);
                fineFluxRegister.incrementFineRegister(fluxes[dd], iter, dd);
                refluxFluxRegister.incrementFineRegister(fluxes[dd], iter, dd);
            }
        }
    }

    template <typename T, unsigned int C, MemType MEM>
    bool computeIncrementFineError(
        MBLevelBoxData<T, C, MEM> &fineRegisters,
        std::map<Point, Box> cfBounds,
        int refRatio,
        Array<T,DIM> gridSpacing,
        int block = -1)
    {
        auto& layout = fineRegisters.layout();
        for (auto iter : layout)
        {
            if (block >= 0 && layout.block(iter) != block) { continue; }
            Box B0 = layout[iter];
            auto &fineData = fineRegisters[iter];
            BoxData<T,C> error(fineData.box());
            error.setVal(0);
            error += fineData;
            // pr_out() << "\nChecking error in patch " << layout.point(iter) << " | " << layout.block(iter) << std::endl;
            // pr_out() << "\tPatch box: " << B0 << std::endl;
            for (auto dir : Point::DirectionsOfCodim(1))
            {
                Box registerBox = cfBounds[dir];
                if (registerBox.empty())
                {
                    continue;
                }
                // pr_out() << "\t\tFound non-empty register in direction " << dir << " | " << registerBox << std::endl;
                auto fineSoln = testFineReflux<T,C>(registerBox, dir, refRatio, gridSpacing);
                // pr_out() << "\t\tComputed solution box: " << fineSoln.box() << std::endl;
                rotateFluxFromFineRegister(fineSoln, iter, dir, layout);
                // pr_out() << "\t\tRotated solution box: " << fineSoln.box() << std::endl;
                error -= fineSoln;
            }
            if (error.absMax() > 1e-12)
            {
                // pr_out() << "\nIncrement Fine Error Detected" << std::endl;
                // pr_out() << "\tblock: " << block << " | patch: " << layout.point(iter) << std::endl;
                // pr_out() << "\tcomputed data: " << std::endl;
                // fineData.printData();
                // pr_out() << "\terror data: " << std::endl;
                // error.printData();
                // BoxData<T,C> trueData(error.box());
                // trueData.setVal(0);
                // trueData += fineData;
                // trueData -= error;
                // pr_out() << "\ttrue data: " << std::endl;
                // trueData.printData();
                return false;
            }
        }
        return true;
    }

    template <typename T, unsigned int C, MemType MEM>
    bool computeIncrementCoarseError(
        MBLevelBoxData<T, C, MEM> &coarseRegisters,
        std::map<Point, Box> cfBounds,
        Array<T,DIM> gridSpacing,
        int block = -1)
    {
        auto& layout = coarseRegisters.layout();
        for (auto iter : layout)
        {
            if (block >= 0 && layout.block(iter) != block) { continue; }
            Box B0 = layout[iter];
            auto &coarseData = coarseRegisters[iter];
            BoxData<T,C> error(coarseData.box());
            error.setVal(0);
            error += coarseData;
            for (auto dir : Point::DirectionsOfCodim(1))
            {
                Box registerBox = cfBounds[dir];
                if (registerBox.empty())
                {
                    continue;
                }
                auto coarseSoln = testCoarseReflux<T,C>(registerBox, dir,gridSpacing);
                error -= coarseSoln;
            }
            if (error.absMax() > 1e-12) { return false; }
        }
        return true;
    }

    template <typename T, unsigned int C, MemType MEM>
    bool computeRefluxError(
        MBLevelBoxData<T, C, MEM> &refluxRegisters,
        std::map<Point, Box> coarseCFBounds,
        int refRatio,
        Array<T,DIM> gridSpacing,
        int block = -1)
    {
        auto& layout = refluxRegisters.layout();
        for (auto iter : layout)
        {
            if (block >= 0 && layout.block(iter) != block) { continue; }
            Box B0 = layout[iter];
            auto &refluxData = refluxRegisters[iter];
            BoxData<T,C> error(refluxData.box());
            error.setVal(0);
            error += refluxData;
            for (auto dir : Point::DirectionsOfCodim(1))
            {
                Box registerBox = coarseCFBounds[dir];
                if (registerBox.empty())
                {
                    continue;
                }
                BoxData<T,C> refluxSoln(registerBox);
                refluxSoln.setVal(0);
                auto coarseSoln = testCoarseReflux<T,C>(registerBox, dir, gridSpacing);
                auto fineArgs = getFineRegisterArgsFromCoarse(dir, registerBox, layout, iter);
                auto fineSoln = testFineReflux<T,C>(get<1>(fineArgs), get<0>(fineArgs), refRatio,gridSpacing);
                if (get<2>(fineArgs) != layout.block(iter))
                {
                    auto R = layout.domain().graph().rotation(layout.block(iter), dir);
                    auto Rinv = R.inverse();
                    fineSoln.rotate(registerBox, Rinv);
                }
                refluxSoln += fineSoln;
                refluxSoln += coarseSoln;

                error -= refluxSoln;
            }
            if (error.absMax() > 1e-12)
            {
                return false;
            }
        }
        return true;
    }


#if DIM == 3
#ifdef PR_AMR
#ifdef PR_MMB
    std::vector<MBPoint> getTestAMRGridPatches(
        Box domainBox,
        Point boxSizes,
        BlockIndex block)
    {
        Box patchDomain = domainBox.coarsen(boxSizes);
        std::vector<MBPoint> patches;
        for (auto patch : patchDomain.edge(Point::X()))
        {
            patches.push_back(MBPoint(patch, block));
        }
        return patches;
    }

    MBAMRGrid testCubedSphereGrid(
        int domainSize,
        int thickness,
        int boxSize, 
        Point refRatio)
    {
        int numBlocks = 6;
        auto domain = CubedSphereShell::Domain(domainSize, thickness, 2);
        Point boxSizeVect(boxSize, boxSize, thickness);
        std::vector<Point> boxSizes(numBlocks, boxSizeVect);
        auto coarsePatches = domain.patches(boxSizeVect);
        MBAMRGrid grid(domain, boxSizeVect, refRatio, 2);

        for (int li = 1; li < 2; li++)
        {
            auto domain = grid[li].domain();
            std::vector<MBPoint> patches;
            for (BlockIndex bi = 0; bi < 6; bi++)
            {
                auto blockPatches = getTestAMRGridPatches(domain.getBlock(bi).box(), boxSizeVect, bi);
                for (auto patch : blockPatches) { patches.push_back(patch); }
            }
            grid[li].define(domain, patches, boxSizes);
        }
        return grid;
    }

    std::map<Point, Box> cubedSphereShellCFBoundary_Coarse(int domainSize, int boxSize, int thickness, int refRatio, BlockIndex block)
    {
        std::map<Point, Box> boundMap;
        Box domainBox(Point(domainSize, domainSize, thickness));
        int coarsenedFinePatchSize = boxSize / refRatio;

        for (auto dir : Point::DirectionsOfCodim(1))
        {
            boundMap[dir] = Box();
        }

        Box unrefinedRegion = domainBox.extrude(Point::X(), -coarsenedFinePatchSize);
        Box posYStrip = domainBox.edge(Point::Y(), coarsenedFinePatchSize);
        Box negYStrip = domainBox.edge(-Point::Y(), coarsenedFinePatchSize);
        Box posXStrip = domainBox.edge(Point::X(), coarsenedFinePatchSize);
        Box negXStrip = domainBox.edge(-Point::X(), coarsenedFinePatchSize);
        boundMap[Point::X()] = unrefinedRegion.edge(0, Side::Hi, 1);

        switch (block)
        {
        case 0:
            boundMap[-Point::X()] = posYStrip.edge(0, Side::Lo, 1);
            boundMap[-Point::Y()] = negXStrip.edge(1, Side::Lo, 1);
            break;
        case 1:
            boundMap[-Point::X()] = negYStrip.edge(0, Side::Lo, 1);
            boundMap[Point::Y()] = negXStrip.edge(1, Side::Hi, 1);
            break;
        case 2:
            boundMap[-Point::X()] = unrefinedRegion.edge(0, Side::Lo, 1);
            break;
        case 3:
            boundMap[-Point::X()] = unrefinedRegion.edge(0, Side::Lo, 1);
            boundMap[Point::Y()] = unrefinedRegion.edge(1, Side::Hi, 1);
            boundMap[-Point::Y()] = unrefinedRegion.edge(1, Side::Lo, 1);
            break;
        case 4:
            boundMap[-Point::X()] = unrefinedRegion.edge(0, Side::Lo, 1);
            boundMap[Point::Y()] = negXStrip.edge(1, Side::Hi, 1);
            boundMap[-Point::Y()] = negXStrip.edge(1, Side::Lo, 1);
            break;
        case 5:
            boundMap[-Point::X()] = unrefinedRegion.edge(0, Side::Lo, 1);
            break;
        }

        return boundMap;
    }
    std::map<Point, Box> cubedSphereShellCFBoundary_Fine(int domainSize, int boxSize, int thickness, int refRatio, BlockIndex block)
    {
        std::map<Point, Box> boundMap;
        Box domainBox(Point(domainSize, domainSize, thickness));
        int coarsenedFinePatchSize = boxSize / refRatio;

        for (auto dir : Point::DirectionsOfCodim(1))
        {
            boundMap[dir] = Box();
        }

        Box refinedRegion = domainBox.edge(Point::X() * coarsenedFinePatchSize);
        boundMap[-Point::X()] = refinedRegion.adjacent(0, Side::Lo, 1);
        boundMap[Point::X()] = refinedRegion.adjacent(0, Side::Hi, 1);

        switch (block)
        {
        case 0:
            boundMap[-Point::Y()] = refinedRegion.adjacent(1, Side::Lo, 1);
            boundMap[Point::X()] = boundMap[Point::X()].extrude(-Point::Y(), -coarsenedFinePatchSize);
            break;
        case 1:
            boundMap[Point::Y()] = refinedRegion.adjacent(1, Side::Hi, 1);
            boundMap[Point::X()] = boundMap[Point::X()].extrude(Point::Y(), -coarsenedFinePatchSize);
            break;
        case 4:
        case 5:
            boundMap[-Point::Y()] = refinedRegion.adjacent(1, Side::Lo, 1);
            boundMap[Point::Y()] = refinedRegion.adjacent(1, Side::Hi, 1);
            break;
        default:
            break;
        }

        return boundMap;
    }

    template <typename T, unsigned int C, MemType MEM>
    bool computeIncrementFineErrorForCubedSphere(
        MBLevelBoxData<T, C, MEM> &fineRegisterData,
        MBLevelFluxRegister<T, C, MEM> &fluxRegister,
        int refRatio,
        Array<T, DIM> gridSpacing)
    {
        auto layout = fineRegisterData.layout();
        MBLevelFluxRegisterTester<T, C, MEM> tester(fluxRegister);
        MBLevelBoxData<T, C, MEM> error(layout, Point::Zeros());
        error.setVal(0);
        for (auto iter : layout)
        {
            auto &fineRegisters = tester.getFineRegistersAtIndex(iter);
            for (auto &ri : fineRegisters)
            {
                auto fineSoln = testFineReflux<T, C>(ri.data().box(), ri.dir(), refRatio, gridSpacing);
                MBIndex adjIndex = rotateFluxFromFineRegister(fineSoln, iter, ri.dir(), layout);
                PROTO_ASSERT(adjIndex != *layout.end(), "Error: Data Corruption");
                PROTO_ASSERT(layout[adjIndex].containsBox(fineSoln.box()), "Error: Data Corruption");

                auto &error_i = error[adjIndex];
                auto &data_i = fineRegisterData[adjIndex];
                BoxData<T, C, MEM> tmpData(fineSoln.box());
                data_i.copyTo(tmpData);

                error_i += tmpData;
                error_i -= fineSoln;
            }
        }
        T errorNorm = error.absMax();
        bool success = (errorNorm < 1e-12);
        return success;
    }

#endif
#endif
#endif
}

#if PR_MMB
#if 1
TEST(MBLevelFluxRegister, TelescopingXPointConstruction) {
    #if PR_VERBOSE > 0
        HDF5Handler h5;
    #endif

    int domainSize = 16;
    int boxSize = 8;
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
            Box cfBound = cfBounds_Coarse[dir] & grid[0][iter];
            bool success = checkRegisters(coarseRegisters, cfBound, dir);
            EXPECT_TRUE(success);
        }
    }

    auto cfBounds_Fine = telescopingCFBoundary_Fine(domainSize);
    for (auto iter : grid[1])
    {
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
    int boxSize = 8;
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
            Box cfBound = cfBounds_Coarse[dir] & grid[0][iter];
            EXPECT_TRUE(checkRegisters(coarseRegisters, cfBound, dir));
        }
    }

    auto cfBounds_Fine = refinedBlockBoundaryCFBoundary_Fine(domainSize, boxSize, refRatio);
    for (auto iter : grid[1])
    {
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

    initializeFluxRegisters(coarseFluxRegister, fineFluxRegister, refluxFluxRegister, grid);

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

    h5.writeMBLevel(map[0], coarseRegisters, "TELESCOPING_COARSE_REGISTERS");
    h5.writeMBLevel(map[0], fineRegisters, "TELESCOPING_FINE_REGISTERS");
    h5.writeMBLevel(map[0], refluxRegisters, "TELESCOPING_REFLUX_REGISTERS");
#endif
    auto coarseCFBounds = telescopingCFBoundary_Coarse(domainSize);
    EXPECT_TRUE(computeIncrementCoarseError(coarseRegisters, coarseCFBounds, gridSpacing));

    auto fineCFBounds = telescopingCFBoundary_Fine(domainSize);
    EXPECT_TRUE(computeIncrementFineError(fineRegisters, fineCFBounds, refRatio, gridSpacing));
    
    EXPECT_TRUE(computeRefluxError(refluxRegisters, coarseCFBounds, refRatio, gridSpacing));
}

TEST(MBLevelFluxRegister, RefinedBlockBoundaryXPointIncrement) {
    constexpr int NUM_COMPS = 1;
    
    int domainSize = 16;
    int boxSize = 16;
    int refRatio = 4;
    int ghostWidth = 2;

    Point refRatios = Point::Ones(refRatio);
    Point ghostWidths = Point::Ones(ghostWidth);

    auto grid = refinedBlockBoundaryXPointGrid(domainSize, 2, refRatio, boxSize);

    Array<double, DIM> gridSpacing = Point::Ones();
    gridSpacing /= domainSize;
    MBLevelFluxRegister<double, NUM_COMPS, HOST> coarseFluxRegister(grid[0], grid[1], refRatios, gridSpacing);
    MBLevelFluxRegister<double, NUM_COMPS, HOST> fineFluxRegister(grid[0], grid[1], refRatios, gridSpacing);
    MBLevelFluxRegister<double, NUM_COMPS, HOST> refluxFluxRegister(grid[0], grid[1], refRatios, gridSpacing);

    initializeFluxRegisters(coarseFluxRegister, fineFluxRegister, refluxFluxRegister, grid);

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

    h5.writeMBLevel(map[0], coarseRegisters, "REFINED_BB_COARSE_REGISTERS");
    h5.writeMBLevel(map[0], fineRegisters, "REFINED_BB_FINE_REGISTERS");
    h5.writeMBLevel(map[0], refluxRegisters, "REFINED_BB_REFLUX_REGISTERS");
#endif
    auto coarseCFBounds = refinedBlockBoundaryCFBoundary_Coarse(domainSize, boxSize, refRatio);
    EXPECT_TRUE(computeIncrementCoarseError(coarseRegisters, coarseCFBounds, gridSpacing));

    auto fineCFBounds = refinedBlockBoundaryCFBoundary_Fine(domainSize, boxSize, refRatio);
    EXPECT_TRUE(computeIncrementFineError(fineRegisters, fineCFBounds, refRatio, gridSpacing));
    
    EXPECT_TRUE(computeRefluxError(refluxRegisters, coarseCFBounds, refRatio, gridSpacing));
}
#endif
#if DIM == 3

TEST(MBLevelFluxRegister, CubedSphereConstruction) {
   
    constexpr int NUMCOMPS = 1;
    int domainSize = 16;
    int boxSize = 8;
    int thickness = 2;
    int ghostWidth = 2;
    int refRatio = 4;

    Point ghostSize(ghostWidth, ghostWidth, 0);
    Point refRatios(refRatio, refRatio, 1);
    Array<double, DIM> gridSpacing;
    gridSpacing.fill(1.0/domainSize);

    MBAMRGrid grid = testCubedSphereGrid(domainSize, thickness, boxSize, refRatios);

    auto& fineLayout = grid[1];
    auto& crseLayout = grid[0];

    auto coarseMap = CubedSphereShell::Map(grid[0], ghostSize);
    auto fineMap = CubedSphereShell::Map(grid[1], ghostSize);
    
    MBLevelFluxRegister<double,NUMCOMPS,HOST> fluxRegister(crseLayout, fineLayout, refRatios, gridSpacing);
    MBLevelFluxRegisterTester<double, NUMCOMPS, HOST> tester(fluxRegister);
    for (BlockIndex bi = 0; bi < 6; bi++)
    {
        auto& localCoarseLayout = crseLayout.getBlock(bi);
        auto cfBoundsCoarse = cubedSphereShellCFBoundary_Coarse(domainSize, boxSize, thickness, refRatio, bi);
        for (auto iter : crseLayout.getBlock(bi))
        {
            
            PatchID patch = localCoarseLayout.point(iter);
            MBIndex globalIter = crseLayout.find(patch, bi);
            PROTO_ASSERT(globalIter != crseLayout.end(),
                "MBLevelFluxRegisterTests | Error: Data Corruption");
            
            auto coarseRegisters = tester.getCoarseRegistersAtIndex(globalIter);
            for (auto dir : Point::DirectionsOfCodim(1))
            {
                Box cfBound = cfBoundsCoarse[dir] & crseLayout[globalIter];
                bool success = checkRegisters(coarseRegisters, cfBound, dir);
                EXPECT_TRUE(success);
            }
        }

        auto &localFineLayout = fineLayout.getBlock(bi);
        auto cfBoundsFine = cubedSphereShellCFBoundary_Fine(domainSize, boxSize, thickness, refRatio, bi);
        for (auto iter : fineLayout.getBlock(bi))
        {

            PatchID patch = localFineLayout.point(iter);
            MBIndex globalIter = fineLayout.find(patch, bi);
            PROTO_ASSERT(globalIter != fineLayout.end(),
                         "MBLevelFluxRegisterTests | Error: Data Corruption");

            auto fineRegisters = tester.getFineRegistersAtIndex(globalIter);
            
            for (auto dir : Point::DirectionsOfCodim(1))
            {
                //pr_out() << "\nchecking fine registers in patch: " << patch << " | dir: " << dir << " | block: " << bi << std::endl;
                Box cfBound = fineLayout[globalIter].adjacent(dir).coarsen(refRatios);
                cfBound &= cfBoundsFine[dir];
                bool success = checkRegisters(fineRegisters, cfBound, dir);
                EXPECT_TRUE(success);
            }
        }
    }

    MBLevelBoxData<double, 1, HOST> fineGrid(grid[1], Point::Zeros());
    MBLevelBoxData<double, 1, HOST> coarseRegisters(grid[0], ghostSize);
    MBLevelBoxData<double, 1, HOST> fineRegisters(grid[0], ghostSize);
    coarseRegisters.setVal(0);
    fineRegisters.setVal(0);
    fineGrid.setVal(0);
    for (BlockIndex bi = 0; bi < 6; bi++)
    {
        auto cfBoundsCoarse = cubedSphereShellCFBoundary_Coarse(domainSize, boxSize, thickness, refRatio, bi);
        auto& levelDataCoarse = coarseRegisters.getBlock(bi);
        for (auto iter : levelDataCoarse)
        {
            auto& patchData = levelDataCoarse[iter];
            for (auto bound : cfBoundsCoarse)
            {
                BoxData<double, 1> boundData(bound.second);
                boundData.setVal(1);
                patchData += boundData;
            }
        }
        auto cfBoundsFine = cubedSphereShellCFBoundary_Fine(domainSize, boxSize, thickness, refRatio, bi);
        auto& levelDataFine = fineRegisters.getBlock(bi);
        for (auto iter : levelDataFine)
        {
            auto& patchData = levelDataFine[iter];
            for (auto bound : cfBoundsFine)
            {
                BoxData<double, 1> boundData(bound.second);
                boundData.setVal(1);
                patchData += boundData;
            }
        }
    }

    HDF5Handler h5;
    h5.writeMBLevel({"data"}, coarseMap, coarseRegisters, "CUBED_SPHERE_COARSE_REGISTERS");
    h5.writeMBLevel({"data"}, coarseMap, fineRegisters, "CUBED_SPHERE_FINE_REGISTERS");
    h5.writeMBLevel({"proxy"}, fineMap, fineGrid, "CUBED_SPHERE_FINE_GRID");

    
}

#if 0 //This test needs reworking
TEST(MBLevelFluxRegister, CubedSphereIncrement) {
    constexpr int NUM_COMPS = 1;
    
    int domainSize = 8;
    int boxSize = 8;
    int refRatio = 4;
    int ghostWidth = 2;
    int thickness = 1;

    Point refRatios = Point::Ones(refRatio);
    Point ghostWidths = Point::Ones(ghostWidth);
    refRatios[2] = 1;
    ghostWidths[2] = 0;

    MBAMRGrid grid = testCubedSphereGrid(domainSize, thickness, boxSize, refRatios);

    Array<double, DIM> gridSpacing = Point::Ones();
    gridSpacing /= domainSize;
    MBLevelFluxRegister<double, NUM_COMPS, HOST> coarseFluxRegister(grid[0], grid[1], refRatios, gridSpacing);
    MBLevelFluxRegister<double, NUM_COMPS, HOST> fineFluxRegister(grid[0], grid[1], refRatios, gridSpacing);
    MBLevelFluxRegister<double, NUM_COMPS, HOST> refluxFluxRegister(grid[0], grid[1], refRatios, gridSpacing);

    initializeFluxRegisters(coarseFluxRegister, fineFluxRegister, refluxFluxRegister, grid);

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

    auto map = CubedSphereShell::Map(grid[0], ghostWidths);

    h5.writeMBLevel(map, coarseRegisters, "REFINED_BB_COARSE_REGISTERS");
    h5.writeMBLevel(map, fineRegisters, "REFINED_BB_FINE_REGISTERS");
    h5.writeMBLevel(map, refluxRegisters, "REFINED_BB_REFLUX_REGISTERS");
#endif

    for (BlockIndex bi = 0; bi < 6; bi++)
    {
        auto coarseCFBounds = cubedSphereShellCFBoundary_Coarse(domainSize, boxSize, thickness, refRatio, bi);
        EXPECT_TRUE(computeIncrementCoarseError(coarseRegisters, coarseCFBounds, gridSpacing, bi));

        //auto fineCFBounds = cubedSphereShellCFBoundary_Fine(domainSize, boxSize, thickness, refRatio, bi);
        //EXPECT_TRUE(computeIncrementFineError(fineRegisters, fineCFBounds, refRatio, gridSpacing, bi));

        //EXPECT_TRUE(computeRefluxError(refluxRegisters, coarseCFBounds, refRatio, gridSpacing, bi));
    }

    EXPECT_TRUE(computeIncrementFineErrorForCubedSphere(fineRegisters, fineFluxRegister, refRatio, gridSpacing));

}
#endif
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
