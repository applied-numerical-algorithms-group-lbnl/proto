#include <gtest/gtest.h>
#include "Proto.H"
#include "TestFunctions.H"
#include "MBMap_Shear.H"

using namespace Proto;


#if DIM == 2
#if 1
TEST(MBLevelMapTests, ShearMap)
{
    int domainSize = 16;
    int boxSize = 8;
    int ghostWidth = 2;
    double gridSpacing = 1.0 / domainSize;
    HDF5Handler h5;

    auto domain = buildShear(domainSize);
    MBDisjointBoxLayout layout(domain, Point::Ones(boxSize));

    // initialize map
    MBLevelMap<MBMap_Shear<HOST>, HOST> map;
    map.define(layout, Point::Ones(ghostWidth));

#if PR_VERBOSE > 0
    h5.writeMBLevel({"X", "Y", "Z"}, map, map.map(), "MBLevelMapTests_ShearMap_X");
#endif
    for (auto iter : layout)
    {
        auto block = layout.block(iter);
        const auto &x_i = map.map()[iter];
        Box box_i = layout[iter].grow(ghostWidth).grow(Centering::PR_NODE);
        EXPECT_EQ(box_i, x_i.box());

        BoxData<double, DIM> x_exact(box_i);
        for (auto point_j : box_i)
        {
            double x_j = point_j[0] * gridSpacing;
            double y_j = point_j[1] * gridSpacing;
            switch (block)
            {
            case 0:
                x_j = x_j - 1.0;
                y_j = y_j - 1.0;
                break;
            case 1:
                y_j = y_j - 1.0 + MB_MAP_SHEAR_SLOPE * x_j;
                break;
            case 2:
                y_j = y_j + MB_MAP_SHEAR_SLOPE * x_j;
                break;
            case 3:
                x_j = x_j - 1.0;
                break;
            }
            x_exact(point_j, 0) = x_j;
            x_exact(point_j, 1) = y_j;
            for (int ii = 2; ii < DIM; ii++)
            {
                x_exact(point_j, ii) = point_j[ii] * gridSpacing;
            }
        }
        auto error_i = x_i - x_exact;

        EXPECT_NEAR(error_i.absMax(), 0, 1e-12);
    }
}
#endif
TEST(MBLevelMapTests, ShearInverseMap)
{
    int domainSize = 16;
    int boxSize = 8;
    int ghostWidth = 2;
    double gridSpacing = 1.0 / domainSize;
    HDF5Handler h5;

    auto domain = buildShear(domainSize);
    MBDisjointBoxLayout layout(domain, Point::Ones(boxSize));

    // initialize map
    MBLevelMap<MBMap_Shear<HOST>, HOST> map;
    map.define(layout, Point::Ones(ghostWidth));

    for (auto iter : layout)
    {
        BlockIndex block = layout.block(iter);
        Box B0 = layout[iter];
        BoxData<double, DIM> X(B0);
        BoxData<double, 1> J(B0);
        map.op(block).apply(X,J);
        BoxData<double, DIM> X0(B0);
        map.op(block).inverse(X0, X);
        auto X0Sln = map.X(B0, gridSpacing);
        X0 -= X0Sln;
        EXPECT_LT(X0.absMax(), 1e-12);
    }
}
TEST(MBLevelMapTests, XPointMapSmall)
{
    int domainSize = 16;
    int boxSize = 8;
    int ghostWidth = 2;
    constexpr int numBlocks = 5;
    double gridSpacing = 1.0 / domainSize;
    HDF5Handler h5;

    auto domain = buildXPoint(domainSize, numBlocks);
    MBDisjointBoxLayout layout(domain, Point::Ones(boxSize));

    // initialize map
    MBLevelMap<MBMap_XPointRigid<numBlocks, HOST>, HOST> map;
    map.define(layout, Point::Ones(ghostWidth));
}
#if 1
TEST(MBLevelMapTests, XPointMap)
{
    int domainSize = 16;
    int boxSize = 8;
    int ghostWidth = 2;
    constexpr int numBlocks = 5;
    double gridSpacing = 1.0 / domainSize;
    HDF5Handler h5;

    auto domain = buildXPoint(domainSize, numBlocks);
    MBDisjointBoxLayout layout(domain, Point::Ones(boxSize));

    // initialize map
    MBLevelMap<MBMap_XPointRigid<numBlocks, HOST>, HOST> map;
    map.define(layout, Point::Ones(ghostWidth));

#if PR_VERBOSE > 0
    h5.writeMBLevel({"X", "Y", "Z"}, map, map.map(), "MBLevelMapTests_XPointMap_X");
#endif
    MBLevelBoxData<double, DIM, HOST> levelError(layout, Point::Zeros());
    MBLevelBoxData<double, DIM, HOST, PR_NODE> levelXExact(layout, Point::Ones(ghostWidth));
    for (auto iter : layout)
    {
        auto block = layout.block(iter);
        const auto &x_i = map.map()[iter];
        Box box_i = layout[iter].grow(ghostWidth).grow(Centering::PR_NODE);
        EXPECT_EQ(box_i, x_i.box());

        auto& x_exact = levelXExact[iter];
        for (auto point_j : box_i)
        {
            double xCenter = 1.0 - point_j[0]*gridSpacing;
            double yCenter = 1.0 - point_j[1]*gridSpacing;

            double dAngleMapped = 2.0*M_PI / numBlocks;
            double dAngleUnmapped = M_PI / 2.0;
            
            double radius = sqrt(xCenter*xCenter + yCenter*yCenter + 2.0*xCenter*yCenter*cos(dAngleMapped));

            double angleAbsolute = M_PI 
                + atan2(yCenter * sin(dAngleMapped), xCenter + yCenter * cos(dAngleMapped)) 
                + dAngleMapped * block;
            
            x_exact(point_j, 0) = radius*cos(angleAbsolute);
            x_exact(point_j, 1) = radius*sin(angleAbsolute);
            for (int ii = 2; ii < DIM; ii++)
            {
                x_exact(point_j, ii) = point_j[ii]*gridSpacing;
            }
        }
        auto& error_i = levelError[iter];
        error_i.setVal(0);
        error_i += x_i;
        error_i -= x_exact;

        EXPECT_NEAR(error_i.absMax(), 0, 1e-12);
    }
    #if PR_VERBOSE > 0
        h5.writeMBLevel({"EX", "EY", "EZ"}, map, levelError, "MBLevelMapTests_XPointMap_Error");
        h5.writeMBLevel({"X", "Y", "Z"}, map, levelXExact, "MBLevelMapTests_XPointMap_XExact");
    #endif
}
TEST(MBLevelMapTests, XPointInverseMap)
{
    int domainSize = 16;
    int boxSize = 16;
    int ghostWidth = 2;
    constexpr int numBlocks = 5;
    double gridSpacing = 1.0 / domainSize;
    HDF5Handler h5;

    auto domain = buildXPoint(domainSize, numBlocks);
    MBDisjointBoxLayout layout(domain, Point::Ones(boxSize));

    // initialize map
    MBLevelMap<MBMap_XPointRigid<numBlocks, HOST>, HOST> map;
    map.define(layout, Point::Ones(ghostWidth));

    for (auto iter : layout)
    {
        BlockIndex block = layout.block(iter);
        Box B0 = layout[iter];
        BoxData<double, DIM> X(B0);
        BoxData<double, 1> J(B0);
        map.op(block).apply(X,J);
        BoxData<double, DIM> X0(B0);
        map.op(block).inverse(X0, X);
        auto X0Sln = map.X(B0, gridSpacing);
        X0 -= X0Sln;
        EXPECT_LT(X0.absMax(), 1e-12);
    }
}
TEST(MBLevelMapTests, InterBlockApply_Shear)
{
    int domainSize = 8;
    int boxSize = 8;
    HDF5Handler h5;

    auto domain = buildShear(domainSize);
    Point boxSizeVect = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);

    Array<Point, DIM + 1> ghost;
    ghost.fill(Point::Ones(4));
    ghost[0] = Point::Ones(1);

    // initialize map
    MBLevelMap<MBMap_Shear<HOST>, HOST> map;
    map.define(layout, ghost);

    for (unsigned int bi = 0; bi < 4; bi++)
    {
        unsigned int bj = (bi + 1) % 4;
        Point dir_ij = layout.domain().graph().connectivity(bi, bj);
        Point dir_ji = layout.domain().graph().connectivity(bj, bi);
        auto Rij = layout.domain().graph().rotation(bi, dir_ij, bj);
        auto Rji = layout.domain().graph().rotation(bj, dir_ji, bi);

        Box Bi = Box::Cube(domainSize).edge(dir_ij, 2);
        Box Bj = Box::Cube(domainSize).adjacent(dir_ji, 2);

        BoxData<double, DIM> Xi(Bi.grow(PR_NODE));
        BoxData<double, DIM> EXi(Bi.grow(PR_NODE));
        BoxData<double, 1> Ji(Bi.grow(PR_NODE));
        BoxData<double, DIM> Xij(Bi.grow(PR_NODE));
        BoxData<double, 1> Jij(Bi.grow(PR_NODE));
        BoxData<double, DIM> Xj(Bj.grow(PR_NODE));
        BoxData<double, DIM> EXj(Bj.grow(PR_NODE));
        BoxData<double, 1> Jj(Bj.grow(PR_NODE));
        BoxData<double, DIM> Xji(Bj.grow(PR_NODE));
        BoxData<double, 1> Jji(Bj.grow(PR_NODE));

        map.apply(Xi, Ji, bi);         // normal compution in bi
        map.apply(Xj, Jj, bj);         // normal compution in bi
        map.doApply(Xij, Jij, bj, bi); // bj func in bi domain
        map.doApply(Xji, Jji, bi, bj); // bi func in bj domain

        Xji.copyTo(EXi, Rji);
        Xij.copyTo(EXj, Rij);

        EXi -= Xi;
        EXj -= Xj;

        EXPECT_LT(EXi.absMax(), 1e-10);
        EXPECT_LT(EXj.absMax(), 1e-10);
    }
}
TEST(MBLevelMapTests, CellApply_Shear)
{
    int domainSize = 8;
    int boxSize = 8;
    HDF5Handler h5;

    auto domain = buildShear(domainSize);
    Point boxSizeVect = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);

    Array<Point, DIM + 1> ghost;
    ghost.fill(Point::Ones(4));
    ghost[0] = Point::Ones(1);

    // initialize map
    MBLevelMap<MBMap_Shear<HOST>, HOST> map;
    map.define(layout, ghost);

    auto XAvg = 0.5 * Shift::Zeros() + 0.5 * Shift::Basis(0);
    auto YAvg0 = 0.5 * Shift::Zeros() + 0.5 * Shift::Basis(1);
    auto YAvg1 = 1.0 * Shift::Basis(1);
    for (auto iter : layout)
    {
        auto block = layout.block(iter);
        BoxData<double, DIM> X_cnr(layout[iter].grow(PR_NODE));
        BoxData<double, 1> J(layout[iter]);
        map.apply(X_cnr, J, block);
        auto X_cnr_0 = slice(X_cnr, 0);
        auto X_cnr_1 = slice(X_cnr, 1);
        BoxData<double, 1> X_ctr_0 = XAvg(X_cnr_0);
        BoxData<double, 1> X_ctr_1(layout[iter]);
        switch (block)
        {
        case 0:
        case 3:
            X_ctr_1 |= YAvg0(X_cnr_1);
            break;
        case 1:
        case 2:
            X_ctr_1 |= YAvg1(X_cnr_1);
            break;
        }
        auto Y_ctr = map.cellCentered(layout[iter], block, block);
        auto Y_ctr_0 = slice(Y_ctr, 0);
        auto Y_ctr_1 = slice(Y_ctr, 1);

        EXPECT_EQ(Y_ctr.box(), layout[iter]);
        Y_ctr_0 -= X_ctr_0;
        Y_ctr_1 -= X_ctr_1;
        EXPECT_LT(Y_ctr_0.absMax(), 1e-12);
        EXPECT_LT(Y_ctr_1.absMax(), 1e-12);
        for (auto pi : Y_ctr.box())
        {
            Array<double, DIM> X_ctr_i;
            X_ctr_i[0] = X_ctr_0(pi, 0);
            X_ctr_i[1] = X_ctr_1(pi, 0);

            MBDataPoint dataPoint(iter, pi, layout);
            auto Y_ctr_i = map.cellCentered(dataPoint);
            Y_ctr_i -= X_ctr_i;
            double err_i = max(abs(Y_ctr_i[0]), abs(Y_ctr_i[1])); // only checking first two dims
            EXPECT_LT(err_i, 1e-12);
        }
    }
}
TEST(MBLevelMapTests, CellApplyBoundary_Shear)
{
    int domainSize = 8;
    int boxSize = 8;
    HDF5Handler h5;

    auto domain = buildShear(domainSize);
    Point boxSizeVect = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);

    Array<Point, DIM + 1> ghost;
    ghost.fill(Point::Ones(4));
    ghost[0] = Point::Ones(1);

    // initialize map
    MBLevelMap<MBMap_Shear<HOST>, HOST> map;
    map.define(layout, ghost);

    MBLevelBoxData<double, 1, HOST> data(layout, ghost);

    for (auto iter : layout)
    {
        auto locBlock = layout.block(iter);
        for (auto dir : Box::Kernel(1))
        {
            auto bounds = data.bounds(iter, dir);
            for (auto bound : bounds)
            {
                Box boundBox = layout[iter].adjacent(ghost[0] * dir);
                auto adjBlock = layout.block(bound.adjIndex);

                auto XLoc = map.cellCentered(boundBox, locBlock, locBlock);
                auto XAdj = map.cellCentered(boundBox, adjBlock, locBlock);
                for (auto bi : boundBox)
                {
                    MBDataPoint bi_loc(iter, bi, layout);
                    MBDataPoint bi_adj(iter, bi, layout, dir, adjBlock);

                    auto XLoc_i = map.cellCentered(bi_loc);
                    auto XAdj_i = map.cellCentered(bi_adj);

                    XLoc_i -= XLoc.array(bi);
                    XAdj_i -= XAdj.array(bi);

                    EXPECT_LT(XLoc_i.absMax(), 1e-12);
                    EXPECT_LT(XAdj_i.absMax(), 1e-12);
                }
            }
        }
    }
}
#endif
#endif
#if DIM > 2
TEST(MBLevelMapTests, CubeSphereShell)
{
    int domainSize = 32;
    int boxSize = 16;
    int thickness = 32;
    int ghostSize = 7;
    int radialDir = CUBED_SPHERE_SHELL_RADIAL_COORD;
    HDF5Handler h5;

    // auto domain = buildCubeSphereShell(domainSize, thickness, radialDir);
    auto domain = CubedSphereShell::Domain(domainSize, thickness, radialDir);
    Point boxSizeVect = Point::Ones(boxSize);
    boxSizeVect[radialDir] = domainSize;
    boxSizeVect[radialDir] = thickness / 2;
    MBDisjointBoxLayout layout(domain, boxSizeVect);
    Point ghost = Point::Ones(ghostSize);

    // initialize map
    MBLevelMap<MBMap_CubedSphereShell<HOST>, HOST> map;
    map.define(layout, ghost);

#if PR_VERBOSE > 0
    h5.writeMBLevel({"X", "Y", "Z"}, map, map.map(), "MBLevelMapTests_CubeSphereMap_X");
#endif
}

TEST(MBLevelMapTests, InterBlockApply_CubeSphereShell)
{
    int domainSize = 8;
    int boxSize = 8;
    int thickness = 8;
    int radialDir = CUBED_SPHERE_SHELL_RADIAL_COORD;
    HDF5Handler h5;

    // auto domain = buildCubeSphereShell(domainSize, thickness, radialDir);
    auto domain = CubedSphereShell::Domain(domainSize, thickness, radialDir);
    Point boxSizeVect = Point::Ones(boxSize);
    boxSizeVect[radialDir] = thickness;
    MBDisjointBoxLayout layout(domain, boxSizeVect);

    Array<Point, DIM + 1> ghost;
    ghost.fill(Point::Ones(4));
    ghost[0] = Point::Ones(0);

    // initialize map
    MBLevelMap<MBMap_CubedSphereShell<HOST>, HOST> map;
    map.define(layout, ghost);

    const auto &graph = layout.domain().graph();
    for (unsigned int bi = 0; bi < 6; bi++)
    {
        for (auto dir_ij : Box::Kernel(1))
        {
            for (auto arc : graph.boundaries(bi, dir_ij))
            {
                unsigned int bj = arc.dstBlock;
                Point dir_ji = arc.dstToSrc;
                auto Rij = arc.R;
                auto Rji = Rij.inverse();

                Box Bi = Box::Cube(domainSize).edge(dir_ij, 2);
                Box Bj = Box::Cube(domainSize).adjacent(dir_ji, 2);

                pr_out() << "Checking boundary between block " << bi << " and " << bj << std::endl;
                pr_out() << "dir_ij: " << dir_ij << " | dir_ji: " << dir_ji << std::endl;

                Box Bi_node = Bi.grow(PR_NODE);
                Box Bj_node = Bj.grow(PR_NODE);
                BoxData<double, DIM> Xi(Bi_node);
                BoxData<double, DIM> EXi(Bi_node);
                BoxData<double, 1> Ji(Bi_node);
                FluxBoxData<double, DIM> NTi(Bi_node);
                BoxData<double, DIM> Xij(Bi_node);
                BoxData<double, 1> Jij(Bi_node);
                BoxData<double, DIM> Xj(Bj_node);
                BoxData<double, DIM> EXj(Bj_node);
                BoxData<double, 1> Jj(Bj_node);
                FluxBoxData<double, DIM> NTj(Bj_node);
                BoxData<double, DIM> Xji(Bj_node);
                BoxData<double, 1> Jji(Bj_node);

                map.apply(Xi, Ji, NTi, bi);    // normal compution in bi
                map.apply(Xj, Jj, NTj, bj);    // normal compution in bi
                map.doApply(Xij, Jij, bj, bi); // bj func in bi domain
                map.doApply(Xji, Jji, bi, bj); // bi func in bj domain

                Xji.copyTo(EXi, Rji);
                Xij.copyTo(EXj, Rij);

                EXi -= Xi;
                EXj -= Xj;

                EXPECT_LT(EXi.absMax(), 1e-10);
                EXPECT_LT(EXj.absMax(), 1e-10);
            }
        }
    }
}
#endif
int main(int argc, char *argv[])
{
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
