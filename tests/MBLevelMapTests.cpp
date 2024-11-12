#include <gtest/gtest.h>
#include "Proto.H"
#include "TestFunctions.H"
#include "MBMap_Shear.H"
//#include "MBMap_CubeSphereShell.H"

using namespace Proto;

#if DIM == 2
TEST(MBLevelMapTests, ShearMap) {
    int domainSize = 8;
    int boxSize = 8;
    HDF5Handler h5;

    auto domain = buildShear(domainSize);
    Point boxSizeVect = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);

    Array<Point, DIM+1> ghost;
    ghost.fill(Point::Ones(4));
    ghost[0] = Point::Ones(1);

    // initialize map
    MBLevelMap<MBMap_Shear, HOST> map;
    map.define(layout, ghost);
    
#if PR_VERBOSE > 0
    h5.writeMBLevel({"X", "Y", "Z"}, map, map.map(), "MBLevelMapTests_ShearMap_X");
#endif
}

TEST(MBLevelMapTests, InterBlockApply_Shear) {
    int domainSize = 8;
    int boxSize = 8;
    HDF5Handler h5;

    auto domain = buildShear(domainSize);
    Point boxSizeVect = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);

    Array<Point, DIM+1> ghost;
    ghost.fill(Point::Ones(4));
    ghost[0] = Point::Ones(1);

    // initialize map
    MBLevelMap<MBMap_Shear, HOST> map;
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
        BoxData<double, 1>   Ji(Bi.grow(PR_NODE));
        BoxData<double, DIM> Xij(Bi.grow(PR_NODE));
        BoxData<double, 1>   Jij(Bi.grow(PR_NODE));
        BoxData<double, DIM> Xj(Bj.grow(PR_NODE));
        BoxData<double, DIM> EXj(Bj.grow(PR_NODE));
        BoxData<double, 1>   Jj(Bj.grow(PR_NODE));
        BoxData<double, DIM> Xji(Bj.grow(PR_NODE));
        BoxData<double, 1>   Jji(Bj.grow(PR_NODE));

        map.apply(Xi, Ji, bi); // normal compution in bi
        map.apply(Xj, Jj, bj); // normal compution in bi
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

    Array<Point, DIM+1> ghost;
    ghost.fill(Point::Ones(4));
    ghost[0] = Point::Ones(1);

    // initialize map
    MBLevelMap<MBMap_Shear, HOST> map;
    map.define(layout, ghost);

    auto XAvg  = 0.5*Shift::Zeros() + 0.5*Shift::Basis(0);
    auto YAvg0 = 0.5*Shift::Zeros() + 0.5*Shift::Basis(1);
    auto YAvg1 = 1.0*Shift::Basis(1);
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
                X_ctr_1 |= YAvg0(X_cnr_1); break;
            case 1:
            case 2:
                X_ctr_1 |= YAvg1(X_cnr_1); break;
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
            double err_i = max(abs(Y_ctr_i[0]), abs(Y_ctr_i[1])); //only checking first two dims
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

    Array<Point, DIM+1> ghost;
    ghost.fill(Point::Ones(4));
    ghost[0] = Point::Ones(1);

    // initialize map
    MBLevelMap<MBMap_Shear, HOST> map;
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
                Box boundBox = layout[iter].adjacent(ghost[0]*dir);
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
#if DIM > 2
TEST(MBLevelMapTests, CubeSphereShell) {
    int domainSize = 32;
    int boxSize = 16;
    int thickness = 32;
    int ghostSize = 7;
    int radialDir = CUBED_SPHERE_SHELL_RADIAL_COORD;
    HDF5Handler h5;

    //auto domain = buildCubeSphereShell(domainSize, thickness, radialDir);
    auto domain = CubedSphereShell::Domain(domainSize, thickness, radialDir);
    Point boxSizeVect = Point::Ones(boxSize);
    boxSizeVect[radialDir] = domainSize;
    boxSizeVect[radialDir] = thickness/2;
    MBDisjointBoxLayout layout(domain, boxSizeVect);
    Point ghost = Point::Ones(ghostSize);

    
    // initialize map
    MBLevelMap<MBMap_CubedSphereShell, HOST> map;
    map.define(layout, ghost);
    
#if PR_VERBOSE > 0
    h5.writeMBLevel({"X", "Y", "Z"}, map, map.map(), "MBLevelMapTests_CubeSphereMap_X");
#endif
}

TEST(MBLevelMapTests, InterBlockApply_CubeSphereShell) {
    int domainSize = 8;
    int boxSize = 8;
    int thickness = 8;
    int radialDir = CUBED_SPHERE_SHELL_RADIAL_COORD;
    HDF5Handler h5;

    //auto domain = buildCubeSphereShell(domainSize, thickness, radialDir);
    auto domain = CubedSphereShell::Domain(domainSize, thickness, radialDir);
    Point boxSizeVect = Point::Ones(boxSize);
    boxSizeVect[radialDir] = thickness;
    MBDisjointBoxLayout layout(domain, boxSizeVect);

    Array<Point, DIM+1> ghost;
    ghost.fill(Point::Ones(4));
    ghost[0] = Point::Ones(0);

    // initialize map
    MBLevelMap<MBMap_CubedSphereShell, HOST> map;
    map.define(layout, ghost);

    const auto& graph = layout.domain().graph();
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
                BoxData<double, 1>   Ji(Bi_node);
                FluxBoxData<double, DIM> NTi(Bi_node);
                BoxData<double, DIM> Xij(Bi_node);
                BoxData<double, 1>   Jij(Bi_node);
                BoxData<double, DIM> Xj(Bj_node);
                BoxData<double, DIM> EXj(Bj_node);
                BoxData<double, 1>   Jj(Bj_node);
                FluxBoxData<double, DIM> NTj(Bj_node);
                BoxData<double, DIM> Xji(Bj_node);
                BoxData<double, 1>   Jji(Bj_node);
                
                map.apply(Xi, Ji, NTi, bi); // normal compution in bi
                map.apply(Xj, Jj, NTj, bj); // normal compution in bi
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
