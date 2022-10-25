#include <gtest/gtest.h>
#include "Proto.H"

#define XPOINT_SIZE 5
#define XPOINT_RADIUS 2.0
#define CUBED_SPHERE_R0 1.0
#define CUBED_SPHERE_R1 2.0

#include "Lambdas.H"

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
bool testCubeSphere(MBLevelBoxData<double, 3, MEMTYPE_DEFAULT, PR_NODE>& a_map, Point a_domainSize, double a_r0, double a_r1)
{
    bool success = true;

    auto& layout = a_map.layout();
    
    Box B0(a_domainSize + Point::Ones());
    Point x = Point::Basis(0);
    Point y = Point::Basis(1);
    Point z = Point::Basis(2);
    
    for (auto iter_0 : layout)
    {
        auto block_0 = layout.block(iter_0);
        auto& patch_0 = a_map[iter_0];
        
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
        auto& D0 = a_map.blockData(b0);
        int b1 = b0+1;
        if (b1 == 6) {b1 = 2;}
        auto& L1 = layout.blockLayout(b1);
        auto& D1 = a_map.blockData(b1);
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
    Array<Point, DIM+1> ghost;
    ghost.fill(Point::Ones());
   
    // requires C++17
    MBMap map(XPointMap, layout, ghost);
    //auto map = buildMap(XPointMap, layout, ghost);

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
    Array<Point, DIM+1> ghost;
    ghost.fill(Point::Ones());
   
    MBMap map(CubedSphereMap, layout, ghost);

    h5.writeMBLevel({"x", "y", "z"}, map.map(), "CUBE_SPHERE.map");
    h5.writeMBLevel({"x", "y", "z"}, map.map(), "CUBE_SPHERE");

    EXPECT_TRUE(testCubeSphere(map.map(), Point::Ones(domainSize), 1.0, 2.0));
}
#endif

#if 0
#define PROTO_MAP_START(name) \
    ACCEL_DECORATION \
    void proto_map_##name (double& x, double& y, double& z, \
            const double& X, const double& Y, const double& Z, unsigned int block) {

#define PROTO_MAP_END(name) \
    } \
    struct struct_proto_map_##name { \
    template<MemType MEM> \
    inline ACCEL_DECORATION \
    void operator() (Point& a_pt, Var<double, 3, MEM>& a_x, \
            unsigned int a_block, \
            Point a_blockSize) const \
    { \
        double X = (a_pt[0]*1.0)/(a_blockSize[0]*1.0); \
        double Y = (a_pt[1]*1.0)/(a_blockSize[1]*1.0); \
        double Z = (a_pt[2]*1.0)/(a_blockSize[2]*1.0); \
        proto_map_##name (a_x(0),a_x(1),a_x(2),X,Y,Z,a_block); \
    } \
    const char* myname = #name; \
    inline Array<double, 3> map(const Array<double, 3>& X, unsigned int block) const { \
        Array<double, 3> x; \
        proto_map_##name (x[0], x[1], x[2], X[0], X[1], X[2], block); \
        return x; \
    } \
    }; \
static struct_proto_map_##name name;

PROTO_MAP_START(mapFoo)
    x = X; y = 2*Y; z = 3*Z;
PROTO_MAP_END(mapFoo)

ACCEL_DECORATION
void f_baz(double& a_output)
{
    a_output = 7;
}

PROTO_KERNEL_START
void f_fooF(Point& a_p, Var<double, 3>& a_data, unsigned int a_block, Point a_blockSize)
{
    f_baz(a_data(0));
}
PROTO_KERNEL_END(f_fooF, f_foo);

template<typename Func>
class TestMap
{
    public:

    TestMap(){};

    TestMap(const Func& a_func)
    {
        m_func = &a_func;
        m_data.define(Box::Cube(9));
        m_data.setVal(0);
    }

    TestMap(const TestMap<Func>& a_map)
    {
        m_func = a_map.m_func;
        m_data.define(Box::Cube(9));
        m_data.setVal(0);
    }
    
    void apply()
    {
        unsigned int block = 0;
        Point blockSize = Point::Ones(8);
        forallInPlace_p(*m_func, m_data, block, blockSize);
    }
    
    Array<double, 3> apply(const Array<double, 3>& a_X, unsigned int block)
    {
        return m_func->map(a_X, block);
    }

    BoxData<double, 3>& data() {return m_data;}
    
    private:
    
    const Func* m_func;
    BoxData<double, 3> m_data;
};

template<typename Func>
TestMap<Func> getMap(const Func& a_func)
{
    TestMap<Func> map(a_func);
    return map;
}

TEST(MBMap, Forall) {
    
    TestMap map(mapFoo);
    map.data().printData();
    map.apply();
    map.data().printData();
    Array<double, 3> X;
    for (int ii = 0; ii < 3; ii++)
    {
        X[ii] = 1.0;
    }
    auto x = map.apply(X, 0);
    for (int ii = 0; ii < 3; ii++)
    {
        std::cout << x[ii] << ", ";
    }
    std::cout << std::endl;
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
