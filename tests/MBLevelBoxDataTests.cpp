#include <gtest/gtest.h>
#include "Proto.H"
#include "TestFunctions.H"

#define NCOMP 1
using namespace Proto;
#if 1
TEST(MBLevelBoxData, Construction) {
    int domainSize = 32;
    int boxSize = 16;
    int numBlocks = 5;
    auto domain = buildXPoint(domainSize, numBlocks);
    Point boxSizeVect = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);
    Array<Point, DIM+1> ghost;
    for (int ii = 0; ii <= DIM; ii++)
    {
        ghost[ii] = Point::Ones(ii+1) + Point::X();
    }
    Point boundGhost = Point::Ones();
    MBLevelBoxData<int, NCOMP, HOST> hostData(layout, ghost, boundGhost);

    Point nx = Point::Basis(0);
    Point ny = Point::Basis(1);
    Box K = Box::Kernel(1);
    Box patchDomain = Box::Cube(domainSize).coarsen(boxSize);
    
    std::map<Point, Point> patchMap;
    for (auto pi : patchDomain.adjacent(nx,1))
    {
        Point p = pi;
        p[0] = pi[1];
        p[1] = patchDomain.size(1)-1;
        patchMap[pi] = p;
    }
    for (auto pi : patchDomain.adjacent(ny,1))
    {
        Point p = pi;
        p[1] = pi[0];
        p[0] = patchDomain.size(0)-1;
        patchMap[pi] = p;
    }
    for (auto pi : patchDomain.adjacent(nx+ny,1))
    {
        patchMap[pi] = pi-(nx+ny);
    }

    auto CW = CoordPermutation::cw();
    auto CCW = CoordPermutation::ccw();

    for (auto iter : layout)
    {
        auto patchID = layout.point(iter);
        auto blockID = layout.block(iter);
        unsigned int xBlock = (blockID+1) % numBlocks;
        unsigned int yBlock = (blockID-1+numBlocks) % numBlocks;
        auto blockLayout = layout.getBlock(blockID);
        Box patchBox = (layout)[iter]; 
        
        for (auto dir : K)
        {
            Point neighbor = patchID + dir;
            Point adjPatch = patchMap[neighbor];
            Box adjPatchBox = Box(adjPatch, adjPatch).refine(boxSize);
            auto bounds = hostData.bounds(iter, dir);
#if PR_VERBOSE > 1
            if (bounds.size() > 0)
            {
                pr_out() << "patch: " << patchID << " | block: " << blockID << " | dir: " << dir << std::endl;
                pr_out() << "\tpatchBox: " << patchBox << " | adjPatchBox: " << adjPatchBox << std::endl;
                pr_out() << "\tnumBoundaries: " << bounds.size() << std::endl;
                for (auto bi : bounds)
                {
                    pr_out() << "\t\tlocalBoundBox: " << bi.localData->box() << " | adjBoundBox: " << bi.adjData->box() << std::endl;
                }
            }
#endif
            if (patchDomain.containsPoint(neighbor))
            {
                EXPECT_EQ(bounds.size(), 0);
            } else if (patchDomain.adjacent(nx,1).containsPoint(neighbor))
            {
                int ghostSize = ghost[1].max();
                EXPECT_EQ(bounds.size(), 1);
                EXPECT_TRUE(hostData.isBlockBoundary(iter, dir, xBlock));
                Box patchBoundary = patchBox.adjacent(dir, ghostSize);
                Point adjDir = -CCW(dir);
                Box adjPatchBoundary = adjPatchBox.edge(adjDir, ghostSize);

                EXPECT_EQ(layout.block(bounds[0].localIndex), blockID);
                EXPECT_EQ(layout.block(bounds[0].adjIndex), xBlock);
                EXPECT_EQ(bounds[0].localData->box(), patchBoundary.grow(boundGhost));
                EXPECT_EQ(bounds[0].adjData->box(), adjPatchBoundary.grow(boundGhost));
            } else if (patchDomain.adjacent(ny,1).containsPoint(neighbor))
            {
                int ghostSize = ghost[1].max();
                EXPECT_EQ(bounds.size(), 1);
                EXPECT_TRUE(hostData.isBlockBoundary(iter, dir, yBlock));
                Box patchBoundary = patchBox.adjacent(dir, ghostSize);
                Point adjDir = -CW(dir); 
                Box adjPatchBoundary = adjPatchBox.edge(adjDir, ghostSize);
                EXPECT_EQ(layout.block(bounds[0].localIndex), blockID);
                EXPECT_EQ(layout.block(bounds[0].adjIndex), yBlock);
                EXPECT_EQ(bounds[0].localData->box(), patchBoundary.grow(boundGhost));
                EXPECT_EQ(bounds[0].adjData->box(), adjPatchBoundary.grow(boundGhost));
            } else if (patchDomain.adjacent(nx+ny,1).containsPoint(neighbor))
            {
                int ghostSize = ghost[2].max();
                EXPECT_EQ(bounds.size(), numBlocks-3);
                Box patchBoundary = patchBox.adjacent(dir,ghostSize);
                Point adjDir = -dir;
                adjDir[0] = dir[0]; adjDir[1] = dir[1];
                Box adjPatchBoundary = adjPatchBox.edge(adjDir, ghostSize);
                for (auto bound : bounds)
                {
                    EXPECT_EQ(layout.block(bound.localIndex), blockID);
                    EXPECT_NE(layout.block(bound.adjIndex), blockID);
                    EXPECT_NE(layout.block(bound.adjIndex), yBlock);
                    EXPECT_NE(layout.block(bound.adjIndex), xBlock);
                    EXPECT_EQ(bound.localData->box(), patchBoundary.grow(boundGhost));
                    EXPECT_EQ(bound.adjData->box(), adjPatchBoundary.grow(boundGhost));
                }
            } else {
                EXPECT_EQ(bounds.size(), 0);
            }
        }
    }
}
#endif
#if 0
TEST(MBLevelBoxData, RefinedConstruction)
{
    int domainSize = 32;
    int boxSize = 16;
    int numBlocks = 5;
    auto domain = buildXPoint(domainSize, numBlocks);
    std::vector<Point> boxSizes;
    std::vector<MBPoint> patches;
    for (BlockIndex block = 0; block < numBlocks; block++)
    {
        boxSizes.push_back(Point::Ones(boxSize));
        patches.push_back(MBPoint(Point::Ones(domainSize / boxSize - 1), block));
    }
    MBDisjointBoxLayout layout(domain, patches, boxSizes);
    MBLevelBoxData<double, 1, HOST> data(layout, Point::Ones(2));
    // TODO: Finish this test
}
#endif
#if 1
TEST(MBLevelBoxData, Initialization) {
    int domainSize = 64;
    int boxSize = 16;
    int ghostSize = 1;
    int numBlocks = 5;
    auto domain = buildXPoint(domainSize, numBlocks);
    Point boxSizeVect = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);

    MBLevelBoxData<double, NCOMP, HOST> hostData(layout, Point::Ones(ghostSize));
    hostData.initialize(f_MBPointID);

    for (auto iter : layout)
    {
        Box patchBox = (layout)[iter];
        auto& hostData_i = hostData[iter];
        int block = layout.block(iter);
        BoxData<double, NCOMP, HOST> tempData(patchBox);
        forallInPlace_p(f_MBPointID, tempData, block);
        for (auto pi : patchBox)
        {
            for (int cc = 0; cc < NCOMP; cc++)
            {
                double err = abs(tempData(pi, cc) - hostData_i(pi, cc));
                EXPECT_LT(err, 1e-12);
            }
        }
    }
}
#endif
#if 1
TEST(MBLevelBoxData, SetBoundary) {
    HDF5Handler h5;
    int domainSize = 64;
    int boxSize = 16;
    int ghostSize = 1;
    int numBlocks = 5;
    auto domain = buildXPoint(domainSize, numBlocks);
    Point boxSizeVect = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);

    MBLevelBoxData<double, DIM, HOST> hostData(layout, Point::Ones(ghostSize));
    hostData.setVal(0);
    for (int ii = 0; ii < DIM; ii++)
    {
        hostData.setBoundary(ii+1,ii);
    }
#if PR_VERBOSE > 0
    h5.writeMBLevel(hostData, "MBLevelBoxDataTests_SetBoundaries"); 
#endif
    for (auto iter : layout)
    {
        auto& patch = hostData[iter];
        for (auto dir : Box::Kernel(1))
        {
            Box b;
            if (dir == Point::Zeros())
            {
                b = (layout)[iter];
            } else {
                b = (layout)[iter].adjacent(dir * ghostSize);
            }
            for (auto pi : b)
            {
                for (int ii = 0; ii < DIM; ii++)
                {
                    if (layout.isDomainBoundary(iter, dir))
                    {
                        EXPECT_EQ(patch(pi, ii), ii+1);
                    } else {
                        EXPECT_EQ(patch(pi, ii), 0);
                    }
                }
            }
        }
    }
    
}
#endif
#if 1
TEST(MBLevelBoxData, Reduce) {
    HDF5Handler h5;
    int domainSize = 16;
    int boxSize = 8;
    int ghostSize = 2;
    int numBlocks = 5;
    auto domain = buildXPoint(domainSize, numBlocks);
    Point boxSizeVect = Point::Ones(boxSize);
    std::vector<Point> boxSizes(numBlocks, boxSizeVect);

    std::vector<MBPoint> patches;
    for (int bi = 0; bi < numBlocks; bi++)
    {
        patches.push_back(MBPoint(Point::Ones(), bi));
    }
    MBDisjointBoxLayout layout(domain, patches, boxSizes);

    Point ghost = Point::Ones(ghostSize);
    ghost += Point::X();
    MBLevelBoxData<double, 1, HOST> hostData(layout, ghost);
    hostData.setVal(-1.0);
    double abs = hostData.absMax();
    double sum = hostData.sum();
    double max = hostData.max();
    double min = hostData.min();

    EXPECT_NEAR(abs, 1, 1e-12);
    double trueSum = 0;
    for (auto iter : layout)
    {
        trueSum += layout[iter].size();
    }
    trueSum *= -1;
    EXPECT_NEAR(sum, trueSum, 1e-12);
    EXPECT_NEAR(max, 0, 1e-12);
    EXPECT_NEAR(min, -1, 1e-12);
}
#endif

#if 1
TEST(MBLevelBoxData, FillBoundaries) {
    HDF5Handler h5;
    int domainSize = 16;
    int boxSize = 8;
    int ghostSize = 2;
    int numBlocks = 5;
    auto domain = buildXPoint(domainSize, numBlocks);
    Point boxSizeVect = Point::Ones(boxSize);
    std::vector<Point> boxSizes(numBlocks, boxSizeVect);
#if 1
    std::vector<MBPoint> patches;
    for (int bi = 0; bi < numBlocks; bi++)
    {
        patches.push_back(MBPoint(Point::Ones(), bi));
    }
    MBDisjointBoxLayout layout(domain, patches, boxSizes);
#else
    MBDisjointBoxLayout layout(domain, boxSizes);
#endif
    Point ghost = Point::Ones(ghostSize);
    ghost += Point::X();
    MBLevelBoxData<double, DIM, HOST> hostData(layout, ghost);
    hostData.setVal(0);
    MBLevelBoxData<double, DIM, HOST> hostData0(layout, Point::Zeros());
    hostData0.initialize(f_MBPointID);
    for (auto iter : layout) { hostData0[iter].copyTo(hostData[iter]); }
#if PR_VERBOSE > 0
    h5.writeMBLevel({"x", "y"}, hostData, "MBLevelBoxData_FillBoundaries_0"); 
#endif
    hostData.exchange();
#if PR_VERBOSE > 0
    h5.writeMBLevel({"x", "y"}, hostData, "MBLevelBoxData_FillBoundaries_1"); 
#endif

    Box dirs = Box::Kernel(1);

    for (auto iter : layout)
    {
        auto block = layout.block(iter);
        auto patch = layout.point(iter);
        for (auto dir : dirs)
        {
            auto bounds = hostData.bounds(iter, dir);
            for (auto bound : bounds)
            {
                auto& localData = *bound.localData;
                BoxData<double, DIM, HOST> adj(bound.adjData->box());
                BoxData<double, DIM, HOST> localSoln(bound.localData->box());
                BoxData<double, DIM, HOST> error(bound.localData->box());
                auto adjBlock = layout.block(bound.adjIndex);
                auto R = bound.adjToLocal;
                forallInPlace_p(f_MBPointID, adj, adjBlock);
                adj.copyTo(localSoln, R);
                localData.copyTo(error);
                error -= localSoln;
                double errNorm = error.absMax();
				EXPECT_LT(errNorm, 1e-12);
#if PR_VERBOSE > 1
				if (errNorm > 1e-12)
				{
					pr_out() << "========================================================" << std::endl;
					pr_out() << "Error in FillBoundaries: block: " << block << " | patch: " << patch << " | dir: " << dir << std::endl;
					pr_out() << "Local Data: " << std::endl;
					localData.printData();
					pr_out() << "Local Solution: " << std::endl;
					localSoln.printData();
					pr_out() << "Error: " << std::endl;
					error.printData();
				}
#endif
			}
        }
    }
}
#endif
#if 1
TEST(MBLevelBoxData, FillBoundaries_Node) {
    HDF5Handler h5;
    int domainSize = 16;
    int boxSize = 8;
    int ghostSize = 2;
    int numBlocks = 5;
    auto domain = buildXPoint(domainSize, numBlocks);
    Point boxSizeVect = Point::Ones(boxSize);
    std::vector<Point> boxSizes(numBlocks, boxSizeVect);

    std::vector<MBPoint> patches;
    for (int bi = 0; bi < numBlocks; bi++)
    {
        patches.push_back(MBPoint(Point::Ones(), bi));
    }
    MBDisjointBoxLayout layout(domain, patches, boxSizes);

    Point ghost = Point::Ones(ghostSize);
    //ghost += Point::X();
    MBLevelBoxData<double, DIM, HOST, PR_NODE> hostData(layout, ghost);
    MBLevelBoxData<double, DIM, HOST, PR_NODE> hostData0(layout, Point::Zeros());
    hostData.setVal(0);
    hostData0.initialize(f_MBPointID);
    for (auto iter : layout) { hostData0[iter].copyTo(hostData[iter]); }
#if PR_VERBOSE > 0
    h5.writeMBLevel({"x", "y"}, hostData, "MBLevelBoxData_FillBoundaries_0"); 
#endif
    hostData.exchange();
#if PR_VERBOSE > 0
    h5.writeMBLevel({"x", "y"}, hostData, "MBLevelBoxData_FillBoundaries_1"); 
#endif

    Box dirs = Box::Kernel(1);

    for (auto iter : layout)
    {
        auto block = layout.block(iter);
        auto patch = layout.point(iter);
        for (auto dir : dirs)
        {
            auto bounds = hostData.bounds(iter, dir);
            for (auto bound : bounds)
            {
                auto& localData = *bound.localData;
                BoxData<double, DIM, HOST> adj(bound.adjData->box());
                BoxData<double, DIM, HOST> localSoln(bound.localData->box());
                BoxData<double, DIM, HOST> error(bound.localData->box());
                auto adjBlock = layout.block(bound.adjIndex);
                auto R = bound.adjToLocal;
                forallInPlace_p(f_MBPointID, adj, adjBlock);
                adj.copyTo(localSoln, R);
                localData.copyTo(error);
                error -= localSoln;
                double errNorm = error.absMax();
				EXPECT_LT(errNorm, 1e-12);
#if PR_VERBOSE > 1
				if (errNorm > 1e-12)
				{
					pr_out() << "========================================================" << std::endl;
					pr_out() << "Error in FillBoundaries: block: " << block << " | patch: " << patch << " | dir: " << dir << std::endl;
					pr_out() << "Local Data: " << std::endl;
					localData.printData();
					pr_out() << "Local Solution: " << std::endl;
					localSoln.printData();
					pr_out() << "Error: " << std::endl;
					error.printData();
				}
#endif
			}
        }
    }
}
#endif
#if 1
TEST(MBLevelBoxData, CopyTo) {
    HDF5Handler h5;
    int domainSize = 128;
    int boxSize = 16;
    int ghostSize = 1;
    int numBlocks = 5;
    auto domain = buildXPoint(domainSize, numBlocks);
    Point boxSizeVect = Point::Ones(boxSize);
    for (int ti = 0; ti < 2; ti++)
    {
        Point srcBoxSize = boxSizeVect;
        Point dstBoxSize = boxSizeVect;
        switch (ti)
        {
            case 0: srcBoxSize *= 2; dstBoxSize *= 4; break;
            case 1: srcBoxSize *= 4; dstBoxSize *= 2; break;
        }
        MBDisjointBoxLayout srcLayout(domain, srcBoxSize);
        MBDisjointBoxLayout dstLayout(domain, dstBoxSize);
        pr_out() << "srcLayout: " << std::endl;
        srcLayout.print();
        pr_out() << "dstLayout: " << std::endl;
        dstLayout.print();
        Point ghost = Point::Ones(ghostSize) + Point::X();
        MBLevelBoxData<double, DIM, HOST> hostSrc(srcLayout, ghost);
        MBLevelBoxData<double, DIM, HOST> hostDst(dstLayout, ghost);
        hostDst.setVal(7);
        hostSrc.initialize(f_MBPointID);
#if PR_VERBOSE > 0
        h5.writeMBLevel({"x", "y", "z"},hostSrc, "CopyTo_Src_T%i_0", ti);
        h5.writeMBLevel({"x", "y", "z"},hostDst, "CopyTo_Dst_T%i_0", ti);
        h5.writePatch({"x", "y", "z"},hostSrc[*srcLayout.begin()], "CopyTo_Src_Patch_T%i_0", ti);
        h5.writePatch({"x", "y", "z"},hostDst[*dstLayout.begin()], "CopyTo_Dst_Patch_T%i_0", ti);
#endif
        hostSrc.copyToFull(hostDst);
#if PR_VERBOSE > 0
        h5.writeMBLevel({"x", "y", "z"},hostSrc, "CopyTo_Src_T%i_1", ti);
        h5.writeMBLevel({"x", "y", "z"},hostDst, "CopyTo_Dst_T%i_1", ti);
        h5.writePatch({"x", "y", "z"},hostSrc[*srcLayout.begin()], "CopyTo_Src_Patch_T%i_1", ti);
        h5.writePatch({"x", "y", "z"},hostDst[*dstLayout.begin()], "CopyTo_Dst_Patch_T%i_1", ti);
#endif

        for (auto iter : dstLayout)
        {
            auto& dst = hostDst[iter];
            auto sln = forall_p<double, DIM>(f_MBPointID, (dstLayout)[iter], dstLayout.block(iter));
            BoxData<double, DIM, HOST> err((dstLayout)[iter]);
            dst.copyTo(err);
            err -= sln;
            EXPECT_LT(err.absMax(), 1e-12);
        }
    }
}
#endif
#if 1
TEST(MBLevelBoxData, OnDomainBoundary)
{
    HDF5Handler h5;
    
    int domainSize = 8;
    int boxSize = 8;
    int numGhost = 2;
    int numBlocks = 5;
    auto domain = buildXPoint(domainSize, numBlocks);
    Point boxSizeVect = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);

    Point ghost = Point::Ones(numGhost);
    MBLevelBoxData<double, DIM, HOST> hostData(layout, ghost);
    hostData.setVal(0);

    for (auto iter : layout)
    {
        auto block = layout.block(iter);
        auto& patch = hostData[iter];

        for (auto pi : patch.box())
        {
            auto domainBoundDirs = hostData.onDomainBoundary(pi, iter);
            for (auto dir : domainBoundDirs)
            {
                for (int di = 0; di < DIM; di++)
                {
                    if (dir[di] != 0)
                    {
                        patch(pi, di) += 1;
                    }
                }
            }
        }
    }
   
    Box validDomain = Box::Cube(domainSize).grow(0,Side::Hi,numGhost).grow(1,Side::Hi, numGhost);
    std::vector<Box> interiorDomains(DIM);
    interiorDomains[0] = validDomain.grow(0,Side::Lo,-1);
    interiorDomains[1] = validDomain.grow(1,Side::Lo,-1);
    for (int di = 2; di < DIM; di++)
    {
        interiorDomains[di] = validDomain.grow(di,-1);
    }
    for (auto iter : layout)
    {
        auto block = layout.block(iter);
        auto& patch = hostData[iter];
        for (auto pi : patch.box())
        {
            for (int dir = 0; dir < DIM; dir++)
            {
                if (!validDomain.containsPoint(pi)) { EXPECT_EQ(patch(pi, dir), 0); }
                else if (interiorDomains[dir].containsPoint(pi)) { EXPECT_EQ(patch(pi, dir), 0); }
                else {
                    EXPECT_EQ(patch(pi,dir), 1);
                }
            }
        }
    }
#if PR_VERBOSE > 0
    std::vector<std::string> varnames;
    for (int dir = 0; dir < DIM; dir++) { varnames.push_back("dir " + std::to_string(dir)); }

    h5.writeMBLevel(varnames, hostData, "OnDomainBoundary_Data");
#endif
}
#endif
#if 1
TEST(MBLevelBoxData, InterpFootprintCorner)
{
    HDF5Handler h5;

    int domainSize = 32;
    int boxSize = 16;
    int numBlocks = 5;
    auto domain = buildXPoint(domainSize, numBlocks);
    Point boxSizeVect = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);

    Array<Point, DIM+1> ghost;
    ghost.fill(Point::Ones(4));
    ghost[0] = Point::Ones(2);

    MBLevelBoxData<double, NCOMP, HOST> hostData(layout, ghost);
    hostData.initialize(f_MBPointID);
    hostData.exchange();

    // input footprint
    std::vector<Point> footprint;
    for (auto pi : Box::Kernel(2))
    {
        if (pi.sumAbs() <= 2)
        {
            footprint.push_back(pi);
        }
    }

    if (procID() == 0)
    {
        // inputs
        Point p0 = Point::Ones(domainSize) + Point::Basis(0);
        Point patchID = Point::Ones((domainSize / boxSize) - 1);
        auto mbIndex = layout.find(patchID, 0);

        // correct output
        Box domainBox_0 = Box::Cube(domainSize);
        Box domainBox_X = domainBox_0.shift(Point::Basis(0,domainSize));
        Box domainBox_Y = domainBox_0.shift(Point::Basis(1,domainSize));
        Box domainBox_XY = domainBox_0.shift(
                Point::Basis(1,domainSize) + Point::Basis(0, domainSize));
        std::vector<MBDataPoint> soln;
        for (auto s : footprint)
        {
            Point p = s + p0;
            if (domainBox_0.containsPoint(p))
            {
                MBDataPoint data(mbIndex, p, layout);
                soln.push_back(data);
            }
            if (domainBox_X.containsPoint(p))
            {
                MBDataPoint data(mbIndex, p, layout, Point::Basis(0), 1);
                soln.push_back(data);
            }
            if (domainBox_Y.containsPoint(p))
            {
                MBDataPoint data(mbIndex, p, layout, Point::Basis(1), numBlocks-1);
                soln.push_back(data);
            }
            if (domainBox_XY.containsPoint(p))
            {
                for (int bi = 2; bi < numBlocks-1; bi++)
                {
                    MBDataPoint data(mbIndex, p, layout, Point::Basis(0) + Point::Basis(1), bi);
                    soln.push_back(data);
                }
            }
        }

        auto mb_footprint = hostData.interpFootprint(p0, ghost[0], footprint, mbIndex);

        EXPECT_EQ(soln.size(), mb_footprint.size());
        for (auto item : soln)
        {
            EXPECT_NE(std::find(mb_footprint.begin(), mb_footprint.end(), item), mb_footprint.end());
        }
    }
}
#endif
#if 1
TEST(MBLevelBoxData, InterpFootprintEdge)
{
    HDF5Handler h5;

    int domainSize = 32;
    int boxSize = 16;
    int numBlocks = 5;
    auto domain = buildXPoint(domainSize, numBlocks);
    Point boxSizeVect = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);

    Array<Point, DIM+1> ghost;
    ghost.fill(Point::Ones(4));
    ghost[0] = Point::Ones(2);

    MBLevelBoxData<double, NCOMP, HOST> hostData(layout, ghost);
    hostData.initialize(f_MBPointID);
    hostData.exchange();

    // input footprint
    std::vector<Point> footprint;
    
    for (auto pi : Box::Kernel(2))
    {
        if (pi.sumAbs() <= 2)
        {
            footprint.push_back(pi);
        }
    }
  
    if (procID() == 0)
    {
        // inputs
        Point p0 = Point::Ones(boxSize/2);
        p0[0] = domainSize;
        p0[1] = boxSize+1;
        //Point p0 = Point::Basis(0,domainSize) + Point::Basis(1,boxSize+1);
        Point patchID = Point::Basis(0, (domainSize / boxSize) - 1);
        auto mbIndex = layout.find(patchID, 0);


        // correct output
        Point nx = Point::Basis(0);
        Point ny = Point::Basis(1);
        Box patchBox_0 = (layout)[mbIndex];
        Box patchBox_X = patchBox_0.adjacent(ghost[1]*nx);
        Box patchBox_XY = patchBox_0.adjacent(ghost[1]*(nx+ny));
        patchBox_0 = patchBox_0.grow(ghost[0]) & Box::Cube(domainSize);
        std::vector<MBDataPoint> soln;
        for (auto s : footprint)
        {
            Point p = s + p0;
            if (patchBox_0.containsPoint(p))
            {
                MBDataPoint data(mbIndex, p, layout);
                soln.push_back(data);
            }
            if (patchBox_X.containsPoint(p))
            {
                MBDataPoint data(mbIndex, p, layout, nx, 1);
                soln.push_back(data);
            }
            if (patchBox_XY.containsPoint(p))
            {
                MBDataPoint data(mbIndex, p, layout, nx + ny, 1);
                soln.push_back(data);
            }
        }

        auto mb_footprint = hostData.interpFootprint(p0, ghost[0], footprint, mbIndex);
        std::sort(mb_footprint.begin(), mb_footprint.end()); 
        std::sort(soln.begin(), soln.end());
        std::sort(footprint.begin(), footprint.end());
#if PR_VERBOSE > 0
        std::cout << "Domain Box: " << Box::Cube(domainSize) << std::endl;
        std::cout << "Destination Point: (" << p0 << " | " << 0 << ") " << std::endl;
        std::cout << "Base Footprint\t\tCorrect Footprint\t\tComputed Footprint" << std::endl;
        int maxLength = max(max(footprint.size(), soln.size()), mb_footprint.size());
        string line = "--------";
        for (int ii = 0; ii < maxLength; ii++)
        {
            if (ii < footprint.size()) { std::cout << footprint[ii] + p0 << " | \t"; } else { std::cout << line; }
            if (ii < soln.size())
            {
                std::cout << " (" << soln[ii].dstPoint() << " | " << soln[ii].block() << ") ";
            } else {
                std::cout << line << line;
            }
            if (ii < mb_footprint.size())
            {
                std::cout << " (" << mb_footprint[ii].dstPoint() << " | " << mb_footprint[ii].block() << ") " << std::endl;
            } else {
                std::cout << line << line << std::endl;
            }
        }
#endif
        EXPECT_EQ(soln.size(), mb_footprint.size());
        for (int ii = 0; ii < soln.size(); ii++)
        {
            EXPECT_EQ(soln[ii].dstPoint(), mb_footprint[ii].dstPoint());
        }
    }
}
#endif
#if DIM == 2
#if 1
TEST(MBLevelBoxData, InterpFootprintDomainBoundary)
{
    HDF5Handler h5;

    int domainSize = 32;
    int boxSize = 16;
    int numBlocks = 5;
    auto domain = buildXPoint(domainSize, numBlocks);
    Point boxSizeVect = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);

    Array<Point, DIM+1> ghost;
    ghost.fill(Point::Ones(4));
    ghost[0] = Point::Ones(2);

    MBLevelBoxData<double, NCOMP, HOST> hostData(layout, ghost);
    hostData.initialize(f_MBPointID);
    hostData.exchange();

    // input footprint
    std::vector<Point> footprint;
    std::vector<Point> footprintExt;
    for (auto pi : Box::Kernel(2))
    {
        if (pi.sumAbs() <= 2)
        {
            footprint.push_back(pi);
            footprintExt.push_back(pi);
        }
    }
    footprintExt.push_back(Point::Y()*3);
    footprintExt.push_back(Point::Y()*2 + Point::X());
    footprintExt.push_back(Point::Y()*2 - Point::X());

    if (procID() == 0)
    {
        // inputs
        Point p0 = Point::Basis(0,domainSize);
        Point patchID = Point::Basis(0,(domainSize / boxSize) - 1);
        auto mbIndex = layout.find(patchID, 0);

        // correct output
        Point nx = Point::Basis(0);
        Point ny = Point::Basis(1);
        Box patchBox_0 = (layout)[mbIndex];
        Box patchBox_X = patchBox_0.adjacent(ghost[1]*nx);
        Box patchBox_XY = patchBox_0.adjacent(ghost[1]*(nx+ny));
        patchBox_0 = patchBox_0.grow(ghost[0]) & Box::Cube(domainSize);
        std::vector<MBDataPoint> soln;
        for (auto s : footprintExt)
        {
            Point p = s + p0;
            if (patchBox_0.containsPoint(p))
            {
                MBDataPoint data(mbIndex, p, layout);
                soln.push_back(data);
            }
            if (patchBox_X.containsPoint(p))
            {
                MBDataPoint data(mbIndex, p, layout, nx, 1);
                soln.push_back(data);
            }
            if (patchBox_XY.containsPoint(p))
            {
                MBDataPoint data(mbIndex, p, layout, nx + ny, 1);
                soln.push_back(data);
            }
        }

        auto mb_footprint = hostData.interpFootprint(p0, ghost[0], footprint, mbIndex);
        std::sort(mb_footprint.begin(), mb_footprint.end()); 
        std::sort(soln.begin(), soln.end());
        EXPECT_EQ(soln.size(), mb_footprint.size());
        for (int ii = 0; ii < soln.size(); ii++)
        {
            EXPECT_EQ(soln[ii].dstPoint(), mb_footprint[ii].dstPoint());
        }
    }
}
#endif
#endif
#if DIM == 3
#if 1
TEST(MBLevelBoxData, InterpFootprintDomainBoundary)
{
    HDF5Handler h5;

    int domainSize = 8;
    int boxSize = 8;
    int thickness = 1;
    int radialDir = CUBED_SPHERE_SHELL_RADIAL_COORD;
    auto domain = CubedSphereShell::Domain(domainSize, thickness, radialDir);
    Point boxSizeVect = Point::Ones(boxSize);
    boxSizeVect[radialDir] = thickness;
    MBDisjointBoxLayout layout(domain, boxSizeVect);

    Array<Point, DIM+1> ghost;
    ghost.fill(Point::Ones(4));
    ghost[0] = Point::Ones(1);

    MBLevelBoxData<double, NCOMP, HOST> hostData(layout, ghost);
    hostData.initialize(f_MBPointID);
    hostData.exchange();

    // input footprint
    std::vector<Point> footprint;
    for (auto pi : Box::Kernel(2))
    {
        if (pi.sumAbs() <= 2)
        {
            footprint.push_back(pi);
        }
    }
    
    // inputs
    if (procID() == 0)
    {
        Point p0 = Point::Basis(0,domainSize);
        Point patchID = Point::Zeros();
        auto mbIndex = layout.find(patchID, 0);

        auto mb_footprint = hostData.interpFootprint(p0, ghost[0], footprint, mbIndex);
    }
}
#endif
#endif
#if 1
TEST(MBLevelBoxData, MBDataPointOperator)
{
    // TODO: Finish implementing this test
    HDF5Handler h5;

    int domainSize = 32;
    int boxSize = 16;
    int numBlocks = 5;
    auto domain = buildXPoint(domainSize, numBlocks);
    Point boxSizeVect = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);

    Array<Point, DIM+1> ghost;
    ghost.fill(Point::Ones(4));
    ghost[0] = Point::Ones(2);

    MBLevelBoxData<double, DIM, HOST> hostData(layout, ghost, Point::Ones());
    hostData.initialize(f_MBPointID);
    hostData.exchange();

    for (auto iter : layout)
    {
        auto block = layout.block(iter);
        for (auto dir : Box::Kernel(1))
        {
            if (dir == Point::Zeros()) { continue; }
            for (auto bound : hostData.bounds(iter, dir))
            {
                auto adjBlock = layout.block(bound.adjIndex);
                for (auto bi : (layout)[iter].adjacent(dir*ghost[0]).grow(1))
                {
                    MBDataPoint pi(iter, bi, layout, dir, adjBlock);
                    for (int ii = 0; ii < DIM; ii++)
                    {
                        //pr_out() << "\tValue at point " << bi << ": " << hostData[pi](ii) << std::endl;
                        double soln = bound.localData->operator()(bi, ii);
                        EXPECT_EQ(soln, hostData[pi](ii));
                    }
                }
            }
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
