#include <gtest/gtest.h>
#include "Proto.H"
#include "Lambdas.H"
using namespace Proto;

DisjointBoxLayout testLayout(int domainSize, Point boxSize)
{
    Box domainBox = Box::Cube(domainSize); 
    Box patchBox = domainBox.coarsen(boxSize);
    std::vector<Point> patches;
    for (auto patch : patchBox)
    {
        if (patch != Point::Zeros()) { patches.push_back(patch); }
    }
    std::array<bool, DIM> periodicity;
    periodicity.fill(true);
    ProblemDomain domain(domainBox, periodicity);
    return DisjointBoxLayout(domain, patches, boxSize);
}

template<typename T, unsigned int C>
bool compareBoxData(
        const BoxData<T, C, HOST>& a_src,
        const BoxData<T, C, HOST>& a_dst)
{
    auto B = a_src.box() & a_dst.box();
    for (auto pt : B)
    {
        for (int cc = 0; cc < C; cc++)
        {
            T diff = abs(a_src(pt, cc) - a_dst(pt, cc));
            if (diff > 1e-12)
            {
                return false;
            }
        }
    }
    return true;
}

template<typename T, unsigned int C>
bool compareLevelData(
        const LevelBoxData<T, C, HOST>& a_src,
        const LevelBoxData<T, C, HOST>& a_dst)
{
    for (auto iter : a_src.layout())
    {
        auto& src = a_src[iter];
        auto& dst = a_dst[iter];
        if (!compareBoxData(src, dst)) { return false; }
    }
    return true; 
}

TEST(LevelBoxData, SetVal) {
    int domainSize = 32;
    double dx = 1.0/domainSize;
    Point boxSize = Point::Ones(16);
    double constVal = 42;
    DisjointBoxLayout layout = testLayout(domainSize, boxSize);
    LevelBoxData<double, 1, HOST> hostData(layout, Point::Ones());
    hostData.setVal(constVal);
#ifdef PROTO_CUDA
    LevelBoxData<double, 1, DEVICE> deviData(layout, Point::Ones());
    deviData.setVal(constVal);
#endif
    for (auto iter : layout)
    {
        auto& hostData_i = hostData[iter];
        int N = hostData_i.size();
        Box B = hostData_i.box();
        BoxData<double, 1, HOST> soln_i(B);
        soln_i.setVal(constVal);
        for (int ii = 0; ii < N; ii++)
        {
            EXPECT_EQ(hostData_i.data()[ii], soln_i.data()[ii]);
        }
#ifdef PROTO_CUDA
        BoxData<double, 1, HOST> tmpData_i(B);
        auto& deviData_i = deviData[iter];
        deviData_i.copyTo(tmpData_i);
        for (int ii = 0; ii < N; ii++)
        {
            EXPECT_EQ(tmpData_i.data()[ii], soln_i.data()[ii]);
        }
#endif
    }
}

TEST(LevelBoxData, Initialize) {
    int domainSize = 32;
    double dx = 1.0/domainSize;
    double offset = 0.125;
    Point boxSize = Point::Ones(16);
    DisjointBoxLayout layout = testLayout(domainSize, boxSize);
    LevelBoxData<double, 1, HOST> hostData(layout, Point::Ones());
    hostData.initialize(f_phi, dx, offset);
#ifdef PROTO_CUDA
    LevelBoxData<double, 1, DEVICE> deviData(layout, Point::Ones());
    deviData.initialize(f_phi, dx, offset);
#endif
    for (auto iter : layout)
    {
        auto& hostData_i = hostData[iter];
        int N = hostData_i.size();
        Box B = hostData_i.box();
        BoxData<double, 1, HOST> soln_i(B);
        forallInPlace_p(f_phi, soln_i, dx, offset);
        EXPECT_TRUE(compareBoxData(soln_i, hostData_i));
#ifdef PROTO_CUDA
        BoxData<double, 1, HOST> tmpData_i(B);
        auto& deviData_i = deviData[iter];
        deviData_i.copyTo(tmpData_i);
        EXPECT_TRUE(compareBoxData(soln_i, tmpData_i));
#endif
    }
}

TEST(LevelBoxData, InitConvolve)
{
    int domainSize = 32;
    double offset = 0.125;
    int numIter = 3;
    Point boxSize = Point::Ones(16);
    double hostErr[numIter];
#ifdef PROTO_CUDA
    double deviErr[numIter];
#endif
    for (int nn = 0; nn < numIter; nn++)
    {
        double dx = 1.0/domainSize;
        auto layout = testLayout(domainSize, boxSize);
        
        LevelBoxData<double, 1, HOST> hostData(layout, Point::Ones());
        LevelBoxData<double, 1, HOST> soln(layout, Point::Ones());
        LevelBoxData<double, 1, HOST> error(layout, Point::Ones());
        hostData.initConvolve(f_phi, dx, offset);
        soln.initialize(f_phi_avg, dx, offset);
        hostErr[nn] = 0;
#ifdef PROTO_CUDA 
        LevelBoxData<double, 1, DEVICE> deviData(layout, Point::Ones());
        deviData.initConvolve(f_phi, dx, offset);
        deviErr[nn] = 0;
#endif
        for (auto iter : layout)
        {
            auto& hostData_i = hostData[iter];
            auto& soln_i = soln[iter];
            auto& error_i = error[iter];
            hostData_i.copyTo(error_i);
            error_i -= soln_i;
        }
        hostErr[nn] = error.absMax();
#ifdef PROTO_CUDA
        for (auto iter : layout)
        {
            auto& deviData_i = deviData[iter];
            auto& soln_i = soln[iter];
            auto& error_i = error[iter];
            deviData_i.copyTo(error_i);
            error_i -= soln_i;
        }
        deviErr[nn] = error.absMax();
#endif
        domainSize *= 2;
    }
    double rate = 4;
    for (int ii = 1; ii < numIter; ii++)
    {
        double hostRate_i = log(hostErr[ii-1]/hostErr[ii])/log(2.0);
        EXPECT_TRUE(abs(rate - hostRate_i) < 0.1);
#ifdef PROTO_CUDA
        double deviRate_i = log(deviErr[ii-1]/deviErr[ii])/log(2.0);
        EXPECT_TRUE(abs(rate - deviRate_i) < 0.1);
#endif
    }
}

TEST(LevelBoxData, LinearSize)
{
    int domainSize = 32;
    double dx = 1.0/domainSize;
    double offset = 0.125;
    double ghostVal = 7;
    Point boxSize = Point::Ones(16);
    auto layout = testLayout(domainSize, boxSize);
    LevelBoxData<double, 1, HOST> src(layout, Point::Ones(1));
    unsigned int srcSize = 0;
    for (auto iter : layout)
    {
        srcSize += src[iter].box().size();
    }
    EXPECT_EQ(srcSize*sizeof(double), src.linearSize());

}
TEST(LevelBoxData, CopyTo)
{
    int domainSize = 64;
    double dx = 1.0/domainSize;
    double offset = 0.125;
    double ghostVal = 7;
    Point boxSize = Point::Ones(16);
    auto layout = testLayout(domainSize, boxSize);
    // HOST -> HOST
    LevelBoxData<double, 2, HOST> hostSrc(layout, Point::Ones(2));
    LevelBoxData<double, 2, HOST> hostDstL(layout, Point::Ones(3));
    LevelBoxData<double, 2, HOST> hostDstS(layout, Point::Ones(1));
    hostSrc.initialize(f_phi, dx, offset);
    hostDstL.setVal(ghostVal);
    hostDstS.setVal(ghostVal);
    hostSrc.copyTo(hostDstL);
    hostSrc.copyTo(hostDstS);
    EXPECT_TRUE(compareLevelData(hostSrc, hostDstL));
    EXPECT_TRUE(compareLevelData(hostSrc, hostDstS));
#ifdef PROTO_CUDA
    // DEVICE -> HOST
    LevelBoxData<double, 2, DEVICE> deviSrc(layout, Point::Ones(2));
    hostDstL.setVal(ghostVal);
    hostDstS.setVal(ghostVal);
    deviSrc.initialize(f_phi, dx, offset);
    hostSrc.initialize(f_phi, dx, offset);
    deviSrc.copyTo(hostDstL);
    deviSrc.copyTo(hostDstS);
    EXPECT_TRUE(compareLevelData(hostSrc, hostDstL));
    EXPECT_TRUE(compareLevelData(hostSrc, hostDstS));
    
    // HOST -> DEVICE
    LevelBoxData<double, 2, DEVICE> deviDstL(layout, Point::Ones(3));
    LevelBoxData<double, 2, DEVICE> deviDstS(layout, Point::Ones(1));
    hostSrc.initialize(f_phi, dx, offset);
    deviDstL.setVal(ghostVal);
    deviDstS.setVal(ghostVal);
    hostSrc.copyTo(deviDstL);
    hostSrc.copyTo(deviDstS);
    for (auto iter : layout)
    {
        auto& deviDstL_i = deviDstL[iter];
        auto& deviDstS_i = deviDstS[iter];
        auto& hostDstL_i = hostDstL[iter];
        auto& hostDstS_i = hostDstS[iter];

        deviDstL_i.copyTo(hostDstL_i);
        deviDstS_i.copyTo(hostDstS_i);
    }
    EXPECT_TRUE(compareLevelData(hostSrc, hostDstL));
    EXPECT_TRUE(compareLevelData(hostSrc, hostDstS));

    // DEVICE -> DEVICE
    hostSrc.initialize(f_phi, dx, offset);
    deviSrc.initialize(f_phi, dx, offset);
    deviDstL.setVal(ghostVal);
    deviDstS.setVal(ghostVal);
    deviSrc.copyTo(deviDstL); 
    deviSrc.copyTo(deviDstS); 
    for (auto iter : layout)
    {
        auto& deviDstL_i = deviDstL[iter];
        auto& deviDstS_i = deviDstS[iter];
        auto& hostDstL_i = hostDstL[iter];
        auto& hostDstS_i = hostDstS[iter];

        deviDstL_i.copyTo(hostDstL_i);
        deviDstS_i.copyTo(hostDstS_i);
    }
    EXPECT_TRUE(compareLevelData(hostSrc, hostDstL));
    EXPECT_TRUE(compareLevelData(hostSrc, hostDstS));
#endif
}

TEST(LevelBoxData, Exchange)
{
    constexpr unsigned int C = 2;
    int domainSize = 64;
    double dx = 1.0/domainSize;
    double offset = 0.125;
    double ghostVal = 7;
    int ghostSize = 1;
    Point boxSize = Point::Ones(16);
    auto layout = testLayout(domainSize, boxSize);
    LevelBoxData<double, C, HOST> hostData(layout, Point::Ones(ghostSize));
    for (auto iter : layout)
    {
        auto& hostData_i = hostData[iter];
        forallInPlace_p(f_phi, hostData_i, dx, offset);
        hostData_i += iter.global();
    }
    hostData.exchange();
    
    for (auto iter : layout)
    {
        auto& hostData_i = hostData[iter];
        Point p = layout.point(iter);
        Box B = layout[iter];
        Box N = Box::Kernel(1).shift(p);
        for (auto n : N)
        {
            if (n == Point::Zeros() || !layout.contains(n)) { continue; }
            Box ghostRegion = B.adjacent(n-p, ghostSize);
            auto ni = layout.find(n);
            BoxData<double, 2, HOST> ghostData(ghostRegion);
            forallInPlace_p(f_phi, ghostData, dx, offset);
            ghostData += ni.global();
            EXPECT_TRUE(compareBoxData(ghostData, hostData_i));
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
