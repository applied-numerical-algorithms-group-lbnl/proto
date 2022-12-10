#include <gtest/gtest.h>
#include "Proto.H"
#include "Lambdas.H"

#define NCOMP 1
using namespace Proto;

TEST(MBLevelBoxData, Construction) {
    int domainSize = 64;
    int boxSize = 16;
    auto domain = buildXPoint(domainSize);
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
        unsigned int xBlock = (blockID+1) % XPOINT_SIZE;
        unsigned int yBlock = (blockID-1+XPOINT_SIZE) % XPOINT_SIZE;
        auto blockLayout = layout.blockLayout(blockID);
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
                EXPECT_EQ(bounds.size(), XPOINT_SIZE-3);
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

TEST(MBLevelBoxData, Initialization) {
    int domainSize = 64;
    int boxSize = 16;
    int ghostSize = 1;
    auto domain = buildXPoint(domainSize);
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

TEST(MBLevelBoxData, FillBoundaries) {
    int domainSize = 2;
    int boxSize = 2;
    int ghostSize = 2;
    auto domain = buildXPoint(domainSize);
    Point boxSizeVect = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);

    MBLevelBoxData<double, NCOMP, HOST> hostData(layout, Point::Ones(ghostSize));
    hostData.initialize(f_MBPointID);
    hostData.fillBoundaries();

    Box dirs = Box::Kernel(1);

    for (auto iter : layout)
    {
        auto block = layout.block(iter);
        auto patch = layout.point(iter);
        pout() << "Block: " << block << std::endl;
        hostData[iter].printData(4);
        hostData.printBounds(iter, 4);
        pout() << "============================================================" << std::endl;
        for (auto dir : dirs)
        {
            auto bounds = hostData.bounds(iter, dir);
            for (auto bound : bounds)
            {
                auto& localData = *bound.localData;
                BoxData<double, NCOMP, HOST> adj(bound.adjData->box());
                BoxData<double, NCOMP, HOST> localSoln(bound.localData->box());
                BoxData<double, NCOMP, HOST> error(bound.localData->box());
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
TEST(MBLevelBoxData, CopyTo) {
    int domainSize = 32;
    int boxSize = 16;
    int ghostSize = 1;
    auto domain = buildXPoint(domainSize);
    Point boxSizeVect = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);

    MBLevelBoxData<double, NCOMP, HOST> hostSrc(layout, Point::Ones(ghostSize));
    MBLevelBoxData<double, NCOMP, HOST> hostDst(layout, Point::Ones(ghostSize));
    hostSrc.initialize(f_MBPointID);
    hostSrc.copyTo(hostDst);

    for (auto iter : layout)
    {
        auto& src = hostSrc[iter];
        auto& dst = hostDst[iter];
        BoxData<double, NCOMP, HOST> err(layout[iter]);
        dst.copyTo(err);
        err -= src;
        EXPECT_LT(err.absMax(), 1e-12);
    }
}

TEST(MBLevelBoxData, InterpFootprintCorner)
{
    HDF5Handler h5;

    int domainSize = 32;
    int boxSize = 16;
    auto domain = buildXPoint(domainSize);
    Point boxSizeVect = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);

    Array<Point, DIM+1> ghost;
    ghost.fill(Point::Ones(4));
    ghost[0] = Point::Ones(2);

    MBLevelBoxData<double, NCOMP, HOST> hostData(layout, ghost);
    hostData.initialize(f_MBPointID);
    hostData.fillBoundaries();

    // input footprint
    std::vector<Point> footprint;
    for (auto pi : Box::Kernel(1))
    {
        footprint.push_back(pi);
    }
    for (int dir = 0; dir < DIM; dir++)
    {
        footprint.push_back(Point::Basis(dir,2));
        footprint.push_back(Point::Basis(dir,-2));
    }

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
            MBDataPoint data(mbIndex, p);
            soln.push_back(data);
        }
        if (domainBox_X.contains(p))
        {
            MBDataPoint data(mbIndex, p, Point::Basis(0), 1);
            soln.push_back(data);
        }
        if (domainBox_Y.contains(p))
        {
            MBDataPoint data(mbIndex, p, Point::Basis(1), XPOINT_SIZE-1);
            soln.push_back(data);
        }
        if (domainBox_XY.contains(p))
        {
            for (int bi = 2; bi < XPOINT_SIZE-1; bi++)
            {
                MBDataPoint data(mbIndex, p, Point::Basis(0) + Point::Basis(1), bi);
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

TEST(MBLevelBoxData, InterpFootprintEdge)
{
    HDF5Handler h5;

    int domainSize = 32;
    int boxSize = 16;
    auto domain = buildXPoint(domainSize);
    Point boxSizeVect = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);

    Array<Point, DIM+1> ghost;
    ghost.fill(Point::Ones(4));
    ghost[0] = Point::Ones(2);

    MBLevelBoxData<double, NCOMP, HOST> hostData(layout, ghost);
    hostData.initialize(f_MBPointID);
    hostData.fillBoundaries();

    // input footprint
    std::vector<Point> footprint;
    for (auto pi : Box::Kernel(1))
    {
        footprint.push_back(pi);
    }
    for (int dir = 0; dir < DIM; dir++)
    {
        footprint.push_back(Point::Basis(dir,2));
        footprint.push_back(Point::Basis(dir,-2));
    }

    // inputs
    Point p0 = Point::Basis(0,domainSize) + Point::Basis(1,boxSize+1);
    Point patchID = Point::Basis(0,(domainSize / boxSize) - 1);
    auto mbIndex = layout.find(patchID, 0);

    // correct output
    Box patchBox_0 = layout[mbIndex];
    Box patchBox_X = patchBox_0.adjacent(Point::Basis(0, ghost[1][0]));
    Box patchBox_XY = patchBox_0.adjacent(ghost[1]);
    patchBox_0 = patchBox_0.grow(ghost[0]) & Box::Cube(domainSize);
    std::vector<MBDataPoint> soln;
    for (auto s : footprint)
    {
        Point p = s + p0;
        if (patchBox_0.contains(p))
        {
            MBDataPoint data(mbIndex, p);
            soln.push_back(data);
        }
        if (patchBox_X.contains(p))
        {
            MBDataPoint data(mbIndex, p, Point::Basis(0), 1);
            soln.push_back(data);
        }
        if (patchBox_XY.contains(p))
        {
            MBDataPoint data(mbIndex, p, Point::Basis(0) + Point::Basis(1), 1);
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
