#include <gtest/gtest.h>
#include "Proto.H"
#include "Lambdas.H"
#include "MBLevelMap_Shear.H"
#include "MBLevelMap_CubeSphereShell.H"

using namespace Proto;

TEST(MBMapTests, ShearMap) {
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
    MBLevelMap_Shear<HOST> map;
    map.define(layout, ghost);

    
#if PR_VERBOSE > 0
    h5.writeMBLevel({"X", "Y", "Z"}, map, map.map(), "MBMapTests_ShearMap_X");
#endif
}

TEST(MBMapTests, InterBlockApply_Shear) {
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
    MBLevelMap_Shear<HOST> map;
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
        BoxData<double, 1>   Ji(Bi);
        BoxData<double, DIM> Xij(Bi.grow(PR_NODE));
        BoxData<double, 1>   Jij(Bi);
        BoxData<double, DIM> Xj(Bj.grow(PR_NODE));
        BoxData<double, DIM> EXj(Bj.grow(PR_NODE));
        BoxData<double, 1>   Jj(Bj);
        BoxData<double, DIM> Xji(Bj.grow(PR_NODE));
        BoxData<double, 1>   Jji(Bj);

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

#if DIM > 2
TEST(MBMapTests, CubeSphereShell) {
    int domainSize = 8;
    int boxSize = 8;
    int thickness = 1;
    int radialDir = 2;
    HDF5Handler h5;

    auto domain = buildCubeSphereShell(domainSize, thickness, radialDir);
    Point boxSizeVect = Point::Ones(boxSize);
    boxSizeVect[radialDir] = thickness;
    MBDisjointBoxLayout layout(domain, boxSizeVect);

    Array<Point, DIM+1> ghost;
    ghost.fill(Point::Ones(4));
    ghost[0] = Point::Ones(1);

    // initialize map
    MBLevelMap_CubeSphereShell<HOST> map;
    map.define(layout, ghost);
    
#if PR_VERBOSE > 0
    h5.writeMBLevel({"X", "Y", "Z"}, map, map.map(), "MBMapTests_CubeSphereMap_X");
#endif
}

TEST(MBMapTests, InterBlockApply_CubeSphereShell) {
    int domainSize = 8;
    int boxSize = 8;
    int thickness = 8;
    int radialDir = 0;
    HDF5Handler h5;

    auto domain = buildCubeSphereShell(domainSize, thickness, radialDir);
    Point boxSizeVect = Point::Ones(boxSize);
    boxSizeVect[radialDir] = thickness;
    MBDisjointBoxLayout layout(domain, boxSizeVect);

    Array<Point, DIM+1> ghost;
    ghost.fill(Point::Ones(4));
    ghost[0] = Point::Ones(0);

    // initialize map
    MBLevelMap_CubeSphereShell<HOST> map;
    map.define(layout, ghost);

    const auto& graph = layout.domain().graph();
    for (unsigned int bi = 0; bi < 1; bi++)
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
                
                pout() << "Checking boundary between block " << bi << " and " << bj << std::endl;
                pout() << "dir_ij: " << dir_ij << " | dir_ji: " << dir_ji << std::endl;

                BoxData<double, DIM> Xi(Bi.grow(PR_NODE));
                BoxData<double, DIM> EXi(Bi.grow(PR_NODE));
                BoxData<double, 1>   Ji(Bi);
                FluxBoxData<double, DIM> NTi(Bi);
                BoxData<double, DIM> Xij(Bi.grow(PR_NODE));
                BoxData<double, 1>   Jij(Bi);
                BoxData<double, DIM> Xj(Bj.grow(PR_NODE));
                BoxData<double, DIM> EXj(Bj.grow(PR_NODE));
                BoxData<double, 1>   Jj(Bj);
                FluxBoxData<double, DIM> NTj(Bj);
                BoxData<double, DIM> Xji(Bj.grow(PR_NODE));
                BoxData<double, 1>   Jji(Bj);
                
                map.apply(Xi, Ji, NTi, bi); // normal compution in bi
                map.apply(Xj, Jj, NTj, bj); // normal compution in bi
                map.doApply(Xij, Jij, bj, bi); // bj func in bi domain
                map.doApply(Xji, Jji, bi, bj); // bi func in bj domain

                Xi.printData();
                Xj.printData();
                Xij.printData();
                Xji.printData();

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
