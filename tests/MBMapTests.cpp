#include <gtest/gtest.h>
#include "Proto.H"
#include "Lambdas.H"

using namespace Proto;


#if DIM > 2
template<typename Func>
bool testCubeSphere(MBMap<Func>& a_map, Point a_domainSize, double a_r0, double a_r1)
{
    bool success = true;

    auto& map = a_map.map();
    auto& layout = map.layout();
    const auto& domain = layout.domain();
    const auto& graph = domain.graph();

    Box B0(a_domainSize + Point::Ones());
    Point x = Point::Basis(0);
    Point y = Point::Basis(1);
    Point z = Point::Basis(2);
    
    for (auto iter_0 : layout)
    {
        auto block_0 = layout.block(iter_0);
        auto& patch_0 = map[iter_0];
        
        Box b_r0 = B0.edge(-z) & patch_0.box();
        Box b_r1 = B0.edge(z)  & patch_0.box();
        for (auto pi : b_r0)
        {
            double r0 = sqrt(patch_0(pi, 0)*patch_0(pi, 0)
                           + patch_0(pi, 1)*patch_0(pi, 1)
                           + patch_0(pi, 2)*patch_0(pi, 2));
            success &= (abs(r0 - a_r0) < 1e-12);
        }
        for (auto pi : b_r1)
        {
            double r1 = sqrt(patch_0(pi, 0)*patch_0(pi, 0)
                           + patch_0(pi, 1)*patch_0(pi, 1)
                           + patch_0(pi, 2)*patch_0(pi, 2));
            success &= (abs(r1 - a_r1) < 1e-12);
        }
    }

    Box Bx_L = B0.edge(-x);
    Box Bx_H = B0.edge(x);
    for (int b0 = 2; b0 < 6; b0++)
    {
        auto& L0 = layout.blockLayout(b0);
        auto& D0 = map.blockData(b0);
        int b1 = b0+1;
        if (b1 == 6) {b1 = 2;}
        auto& L1 = layout.blockLayout(b1);
        auto& D1 = map.blockData(b1);
        for (auto i0 : L0)
        {
            for (auto i1 : L1)
            {
                Point p0 = L0.point(i0);
                Point p1 = L1.point(i1);
                if (!L0.patchDomain().box().edge(x).contains(p0)) {continue; }
                if (!L1.patchDomain().box().edge(-x).contains(p1)) {continue; }
                if (p0[1] != p1[1]) {continue; }
                if (p0[2] != p1[2]) {continue; }
                auto& d0 = D0[i0];
                auto& d1 = D1[i1];
                Box b = B0.edge(-x) & d1.box();
                BoxData<double, 3, HOST> tmp0(b), tmp1(b);
                d0.copyTo(tmp0, d0.box(), -x*a_domainSize);
                d1.copyTo(tmp1);
                tmp0 -= tmp1;
                success &= tmp0.absMax() < 1e-12;
            }
        }
    }

    for (auto iter : layout)
    {
        auto& patch = map[iter];
        auto block = layout.block(iter);
        for (auto pi : patch.box())
        {
            Array<double, 3> x = patch.array(pi);
            Array<double, 3> X;
            for (int xi = 0; xi < 3; xi++)
            {
                X[xi] = (1.0*pi[xi])/(1.0*a_domainSize[xi]);
            }
            auto x_map = a_map(X, block);
            auto x_err = x_map - x;
            double err = x_err.absMax();
            EXPECT_LT(err, 1e-12);
        }

        for (auto pi : patch.box().extrude(Point::Ones(), -1))
        {
            Array<double, 3> x = patch.array(pi + Point::Ones());
            auto x_map = a_map(pi, Point::Ones(), block);
            auto x_err = x_map - x;
            double err = x_err.absMax();
            EXPECT_LT(err, 1e-12);
        }
    }

    for (unsigned int bi = 0; bi < 6; bi++)
    {
        Point dir_ij = Point::Basis(0);
        unsigned int bj = graph.adjacent(bi, dir_ij);
        Box box_0 = domain.blockDomain(bi).box();
        Box box_i = box_0.shift(dir_ij*(a_domainSize / 2));
        Point dir_ji = graph.connectivity(bj, bi);
        Box box_j = box_0.shift(dir_ji*(a_domainSize / 2));

        auto data_i = a_map(box_i, bi, bi);
        auto data_j = a_map(box_i, bi, bj);

        auto soln_i = forall_p<double, 3>(CubedSphereMap, box_i, bi, domain.blockSize(bi), PR_NODE);
        BoxData<double, 3> soln_j(box_j);
        soln_i.copyTo(soln_j, graph.rotation(bi, dir_ij, bj));
       
        BoxData<double, 3> err_i(box_i);
        BoxData<double, 3> err_j(box_j);
        data_i.copyTo(err_i);
        err_i -= soln_i;
        data_j.copyTo(err_j);
        err_j -= soln_j;
       
        double ei = err_i.absMax();
        double ej = err_j.absMax();

        EXPECT_LT(ei, 1e-12);
        EXPECT_LT(ej, 1e-12);
    }
    HDF5Handler h5;
    for (unsigned int bi = 0; bi < 6; bi++)
    {
        Box B = domain.blockDomain(bi).box();
        auto data = a_map(B, bi, PR_CELL);
        double err = 0.0;
        for (auto pi : B)
        {
            Array<double, DIM> x = pi;
            x += 0.5;
            Array<double, DIM> dx = a_domainSize;
            x /= dx;
            auto soln = a_map(x, bi);
            for (int dir = 0; dir < DIM; dir++)
            {
                err = max(err, abs(data(pi, dir) - soln[dir]));
            }
        }
        EXPECT_LT(err, 1e-12);
    }

    auto& J = a_map.jacobian();
    for (auto iter : layout)
    {
        auto& Ji = J[iter];
        EXPECT_GT(Ji.min(), 0);
    }

    return success;
}
#endif

TEST(MBMap, Identity) {
    HDF5Handler h5;
    int domainSize = 8;
    int boxSize = 4;
    auto domain = buildIdentity(domainSize);
    Point boxSizeVect = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);
    Array<Point, DIM+1> ghost;
    ghost.fill(Point::Zeros());
    Point boundGhost = Point::Ones();
   
    // requires C++17
    MBMap<IdentityMap_t> map(IdentityMap, layout, ghost, boundGhost);
    
    auto& J = map.jacobian();
    double dx = 1.0/domainSize;
    double dv = pow(dx,DIM);
    for (auto iter : layout)
    {
        auto& Ji = J[iter];
        // J should be 1/(dx^DIM) everywhere
        EXPECT_LT(abs(Ji.max() - Ji.min()), 1e-12);
        EXPECT_LT(Ji.max()*pow(dx,DIM) - 1, 1e-12);
        EXPECT_GT(Ji.min(), 0);
    }
    for (unsigned int bi = 0; bi < 1; bi++)
    {
        Box B = domain.blockDomain(bi).box();
        auto data = map(B, bi, PR_CELL);
        double err = 0.0;
        for (auto pi : B)
        {
            Array<double, DIM> x = pi;
            x += 0.5;
            x *= dx;
            auto soln = map(x, bi);
            for (int dir = 0; dir < DIM; dir++)
            {
                err = max(err, abs(data(pi, dir) - soln[dir]));
            }
        }
        EXPECT_LT(err, 1e-12);
    }

    h5.writeMBLevel({"x", "y", "z"}, map.map(), "IDENTITY.map");
    h5.writeMBLevel({"J"}, map.jacobian(), "IDENTITY");
}

TEST(MBMap, XPoint) {
    HDF5Handler h5;
    int domainSize = 8;
    int boxSize = 4;
    auto domain = buildXPoint(domainSize);
    Point boxSizeVect = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);
    Array<Point, DIM+1> ghost;
    ghost.fill(Point::Zeros());
    Point boundGhost = Point::Ones();
   
    MBMap<XPointMapRigid_t> map(XPointMapRigid, layout, ghost, boundGhost);
    //auto map = buildMap(XPointMap, layout, ghost);
    auto& J = map.jacobian();
    double dx = 1.0/domainSize;
    double dv = pow(dx,DIM);
    for (auto iter : layout)
    {
        auto& Ji = J[iter];
        EXPECT_GT(Ji.min(), 0);
    }

    h5.writeMBLevel({"x", "y", "z"}, map.map(), "XPOINT.map");
    h5.writeMBLevel({"J"}, map.jacobian(), "XPOINT");
}

#if DIM > 2
TEST(MBMap, CubeSphere) {
    constexpr int C = 2;
    HDF5Handler h5;
    int domainSize = 8;
    int boxSize = 4;
    auto domain = buildCubeSphere(domainSize);
    Point boxSizeVect = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);
    Array<Point, DIM+1> ghost;
    ghost.fill(Point::Zeros());
    Point boundGhost = Point::Ones();
   
    MBMap<CubedSphereMap_t> map(CubedSphereMap, layout, ghost, boundGhost);

    h5.writeMBLevel({"x", "y", "z"}, map.map(), "CUBE_SPHERE.map");
    h5.writeMBLevel({"J"}, map.jacobian(), "CUBE_SPHERE");

    EXPECT_TRUE(testCubeSphere(map, Point::Ones(domainSize), 1.0, 2.0));
}
#endif

TEST(MBMap, InitializeWithMap)
{
    HDF5Handler h5;
    int domainSize = 128;
    int boxSize = 64;
#if DIM > 2
    domainSize = 8;
    boxSize = 4;
#endif
    auto domain = buildXPoint(domainSize);
    Point boxSizeVect = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);
    Array<Point, DIM+1> ghost;
    ghost.fill(Point::Zeros());
    Point boundGhost = Point::Ones();
   
    MBMap<XPointMap_t> map(XPointMap, layout, ghost, boundGhost);

    Point k{1,2,3,4,5,6};
    Array<double, DIM> offset{1,1,1,1,1,1};
    MBLevelBoxData<double, 1, HOST, PR_CELL> hostData(layout, ghost);
    hostData.initialize(f_phiM, map, k, offset);

    h5.writeMBLevel({"x", "y", "z"}, map.map(), "MAP_INIT.map");
    h5.writeMBLevel({"phi"}, hostData, "MAP_INIT");
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
