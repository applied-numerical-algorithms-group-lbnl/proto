#include <gtest/gtest.h>
#include "Proto.H"
#include "Lambdas.H"

#define XPOINT_SIZE 5

using namespace Proto;

#if DIM > 2
MBProblemDomain buildCubeSphere(int a_domainSize)
{
    
    MBProblemDomain domain(6);
    auto R_theta = CoordPermutation::ccw(1);
    auto R_north = CoordPermutation::cw(0);
    auto R_south = CoordPermutation::ccw(0);
    CoordPermutation R_I;
    Point x = Point::Basis(0);
    Point y = Point::Basis(1);
    for (int bi = 2; bi < 6; bi++)
    {
        int srcBlock = bi;
        int dstBlock = bi+1;
        if (dstBlock > 5) { dstBlock = 2; }
        domain.defineBoundary(srcBlock, dstBlock, x, R_I);
        domain.defineBoundary(srcBlock, 1, y, R_north);
        domain.defineBoundary(srcBlock, 0, -y, R_south);
        R_north = R_north*R_theta;
        R_south = R_south*R_theta;
    }
    for (int bi = 0; bi < 6; bi++)
    {
        domain.defineDomain(bi, Point::Ones(a_domainSize));
    }
    return domain;
}
#endif

MBProblemDomain buildXPoint(int a_domainSize)
{
    MBProblemDomain domain(XPOINT_SIZE);
    auto CCW = CoordPermutation::ccw();
    for (int ii = 0; ii < XPOINT_SIZE; ii++)
    {
        domain.defineBoundary(ii, (ii+1) % XPOINT_SIZE, 0, Side::Hi, CCW);
    }
    for (int bi = 0; bi < XPOINT_SIZE; bi++)
    {
        domain.defineDomain(bi, Point::Ones(a_domainSize));
    }
    return domain;
}

#if DIM > 2
bool testCubeSphere(MBMap<double, HOST>& a_map, Point a_domainSize, double a_r0, double a_r1)
{
    bool success = true;

    const auto& map = a_map.map();
    auto& layout = map.layout();
    
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

    return success;
}
#endif

TEST(MBMap, XPoint) {
    pout() << "THIS IS A VISUAL TEST. CONFIRM RESULTS IN VISIT" << std::endl;
    constexpr int C = 2;
    HDF5Handler h5;
    int domainSize = 64;
    int boxSize = 16;
    auto domain = buildXPoint(domainSize);
    Point boxSizeVect = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);
    std::array<Point, DIM+1> ghost;
    ghost.fill(Point::Ones());
    MBLevelBoxData<double, C, HOST> data(layout, Point::Ones());
    data.initialize(f_MBPointID);
    
    MBMap<double, HOST> map(layout, ghost);
    map.compute(f_XPointMap, XPOINT_SIZE, domainSize);
    
    h5.writeMBLevel({"x", "y", "z"}, map.map(), "XPOINT.map");
    h5.writeMBLevel({"x", "y", "z"}, map.map(), "XPOINT");
}

#if DIM > 2
TEST(MBMap, CubeSphere) {
    constexpr int C = 2;
    HDF5Handler h5;
    int domainSize = 64;
    int boxSize = 16;
    auto domain = buildCubeSphere(domainSize);
    Point boxSizeVect = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);
    std::array<Point, DIM+1> ghost;
    ghost.fill(Point::Ones());
    MBLevelBoxData<double, C, HOST> data(layout, Point::Ones());
    data.initialize(f_MBPointID);
    
    MBMap<double, HOST> map(layout, ghost);
    map.compute(f_CubeSphereMap, Point::Ones(domainSize), 1.0, 2.0);
    
    h5.writeMBLevel({"x", "y", "z"}, map.map(), "CUBE_SPHERE.map");
    h5.writeMBLevel({"x", "y", "z"}, data, "CUBE_SPHERE");

    EXPECT_TRUE(testCubeSphere(map, Point::Ones(domainSize), 1.0, 2.0));
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
