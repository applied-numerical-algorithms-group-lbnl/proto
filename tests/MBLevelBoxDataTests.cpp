#include <gtest/gtest.h>
#include "Proto.H"
#include "Lambdas.H"

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
    ghost.fill(Point::Ones());
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
        Box patchBox = layout[iter]; 
        for (auto dir : K)
        {
            Point neighbor = patchID + dir;
            Point adjPatch = patchMap[neighbor];
            Box adjPatchBox = Box(adjPatch, adjPatch).refine(boxSize);
            auto bounds = hostData.bounds(iter, dir);
            if (patchDomain.contains(neighbor))
            {
                EXPECT_EQ(bounds.size(), 0);
            } else if (patchDomain.adjacent(nx,1).contains(neighbor))
            {
                EXPECT_EQ(bounds.size(), 1);
                EXPECT_TRUE(hostData.isBlockBoundary(iter, dir, xBlock));
                Box patchBoundary = patchBox.adjacent(dir, 1);
                Point adjDir = -CCW(dir);
                Box adjPatchBoundary = adjPatchBox.edge(adjDir, 1);

                EXPECT_EQ(layout.block(bounds[0].localIndex), blockID);
                EXPECT_EQ(layout.block(bounds[0].adjIndex), xBlock);
                EXPECT_EQ(bounds[0].localData->box(), patchBoundary.grow(boundGhost));
                EXPECT_EQ(bounds[0].adjData->box(), adjPatchBoundary.grow(boundGhost));
            } else if (patchDomain.adjacent(ny,1).contains(neighbor))
            {
                EXPECT_EQ(bounds.size(), 1);
                EXPECT_TRUE(hostData.isBlockBoundary(iter, dir, yBlock));
                Box patchBoundary = patchBox.adjacent(dir, 1);
                Point adjDir = -CW(dir); 
                Box adjPatchBoundary = adjPatchBox.edge(adjDir, 1);
                EXPECT_EQ(layout.block(bounds[0].localIndex), blockID);
                EXPECT_EQ(layout.block(bounds[0].adjIndex), yBlock);
                EXPECT_EQ(bounds[0].localData->box(), patchBoundary.grow(boundGhost));
                EXPECT_EQ(bounds[0].adjData->box(), adjPatchBoundary.grow(boundGhost));
            } else if (patchDomain.adjacent(nx+ny,1).contains(neighbor))
            {
                EXPECT_EQ(bounds.size(), numBlocks-3);
                Box patchBoundary = patchBox.adjacent(dir,1);
                Point adjDir = -dir;
                adjDir[0] = dir[0]; adjDir[1] = dir[1];
                Box adjPatchBoundary = adjPatchBox.edge(adjDir, 1);
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
        Box patchBox = layout[iter];
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
TEST(MBLevelBoxData, FillBoundaries) {
    int domainSize = 32;
    int boxSize = 32;
    int ghostSize = 4;
    int numBlocks = 5;
    auto domain = buildXPoint(domainSize, numBlocks);
    Point boxSizeVect = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);

    MBLevelBoxData<double, DIM, HOST> hostData(layout, Point::Ones(ghostSize));
    hostData.initialize(f_MBPointID);
    hostData.exchange();

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
            }
        }
    }
}
#endif
TEST(MBLevelBoxData, CopyTo) {
    HDF5Handler h5;
    int domainSize = 32;
    int boxSize = 16;
    int ghostSize = 1;
    int numBlocks = 5;
    auto domain = buildXPoint(domainSize, numBlocks);
    Point boxSizeVect = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);

    MBLevelBoxData<double, DIM, HOST> hostSrc(layout, Point::Ones(ghostSize));
    MBLevelBoxData<double, DIM, HOST> hostDst(layout, Point::Ones(ghostSize));
    hostSrc.initialize(f_MBPointID);
    hostSrc.copyTo(hostDst);

#if PR_VERBOSE > 0
    for (auto iter : layout)
    {
        int block = layout.block(iter);
        h5.writeLevel(1, hostSrc.blockData(block), "CopyTo_Src_B%i", block); 
        h5.writeLevel(1, hostDst.blockData(block), "CopyTo_Dst_B%i", block); 
    }
#endif

    for (auto iter : layout)
    {
        auto& src = hostSrc[iter];
        auto& dst = hostDst[iter];
        BoxData<double, DIM, HOST> err(layout[iter]);
        dst.copyTo(err);
        err -= src;
        EXPECT_LT(err.absMax(), 1e-12);
    }
}
TEST(MBLevelBoxData, OnDomainBoundary)
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
    
    for (auto iter : layout)
    {
        auto block = layout.block(iter);
        auto& patch = hostData[iter];

        for (auto pi : patch.box())
        {
            for (int dir = 0; dir < DIM; dir++)
            {
                if (pi[dir] == 0)
                {
                    EXPECT_EQ(patch(pi,dir), 1);
                } else {
                    EXPECT_EQ(patch(pi,dir), 0);
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
            if (domainBox_0.contains(p))
            {
                MBDataPoint data(mbIndex, p, layout);
                soln.push_back(data);
            }
            if (domainBox_X.contains(p))
            {
                MBDataPoint data(mbIndex, p, layout, Point::Basis(0), 1);
                soln.push_back(data);
            }
            if (domainBox_Y.contains(p))
            {
                MBDataPoint data(mbIndex, p, layout, Point::Basis(1), numBlocks-1);
                soln.push_back(data);
            }
            if (domainBox_XY.contains(p))
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
        Point p0 = Point::Basis(0,domainSize) + Point::Basis(1,boxSize+1);
        Point patchID = Point::Basis(0,(domainSize / boxSize) - 1);
        auto mbIndex = layout.find(patchID, 0);

        // correct output
        Point nx = Point::Basis(0);
        Point ny = Point::Basis(1);
        Box patchBox_0 = layout[mbIndex];
        Box patchBox_X = patchBox_0.adjacent(ghost[1]*nx);
        Box patchBox_XY = patchBox_0.adjacent(ghost[1]*(nx+ny));
        patchBox_0 = patchBox_0.grow(ghost[0]) & Box::Cube(domainSize);
        std::vector<MBDataPoint> soln;
        for (auto s : footprint)
        {
            Point p = s + p0;
            if (patchBox_0.contains(p))
            {
                MBDataPoint data(mbIndex, p, layout);
                soln.push_back(data);
            }
            if (patchBox_X.contains(p))
            {
                MBDataPoint data(mbIndex, p, layout, nx, 1);
                soln.push_back(data);
            }
            if (patchBox_XY.contains(p))
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
            EXPECT_EQ(soln[ii].point, mb_footprint[ii].point);
        }
    }
}
#if DIM == 2
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
        Box patchBox_0 = layout[mbIndex];
        Box patchBox_X = patchBox_0.adjacent(ghost[1]*nx);
        Box patchBox_XY = patchBox_0.adjacent(ghost[1]*(nx+ny));
        patchBox_0 = patchBox_0.grow(ghost[0]) & Box::Cube(domainSize);
        std::vector<MBDataPoint> soln;
        for (auto s : footprintExt)
        {
            Point p = s + p0;
            if (patchBox_0.contains(p))
            {
                MBDataPoint data(mbIndex, p, layout);
                soln.push_back(data);
            }
            if (patchBox_X.contains(p))
            {
                MBDataPoint data(mbIndex, p, layout, nx, 1);
                soln.push_back(data);
            }
            if (patchBox_XY.contains(p))
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
            EXPECT_EQ(soln[ii].point, mb_footprint[ii].point);
        }
    }
}
#endif
#if DIM == 3
TEST(MBLevelBoxData, InterpFootprintDomainBoundary)
{
    HDF5Handler h5;

    int domainSize = 8;
    int boxSize = 8;
    int thickness = 1;
    int radialDir = 2;
    auto domain = buildCubeSphereShell(domainSize, thickness, radialDir);
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
                for (auto bi : layout[iter].adjacent(dir*ghost[0]).grow(1))
                {
                    MBDataPoint pi(iter, bi, layout, dir, adjBlock);
                    for (int ii = 0; ii < DIM; ii++)
                    {
                        //pout() << "\tValue at point " << bi << ": " << hostData[pi](ii) << std::endl;
                        double soln = bound.localData->operator()(bi, ii);
                        EXPECT_EQ(soln, hostData[pi](ii));
                    }
                }
            }
        }
    }
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
