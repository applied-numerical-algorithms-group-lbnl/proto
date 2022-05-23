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
bool testCopyTo(
        const LevelBoxData<T, C, HOST>& a_src,
        const LevelBoxData<T, C, HOST>& a_dst)
{
    for (auto iter : a_src.layout())
    {
        auto& src = a_src[iter];
        auto& dst = a_dst[iter];
        for (auto pt : dst.box())
        {
            for (int cc = 0; cc < C; cc++)
            {
                if (src.box().contains(pt))
                {
                    T diff = abs(src(pt, cc) - dst(pt, cc));
                    if (diff > 1e-12)
                    {
                        return false;
                    }
                }
            }
        }
    }
    return true; 
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
            auto& deviData_i = deviData[lvl][iter];
            auto& soln_i = soln[lvl][iter];
            auto& error_i = error[lvl][iter];
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
    HDF5Handler h5;
    int domainSize = 64;
    double dx = 1.0/domainSize;
    double offset = 0.125;
    double ghostVal = 7;
    Point boxSize = Point::Ones(16);
    auto layout = testLayout(domainSize, boxSize);
    LevelBoxData<double, 1, HOST> hostSrc(layout, Point::Ones(2));
    LevelBoxData<double, 1, HOST> hostDstL(layout, Point::Ones(3));
    LevelBoxData<double, 1, HOST> hostDstS(layout, Point::Ones(1));
    hostSrc.initialize(f_phi, dx, offset);
    hostDstL.setVal(ghostVal);
    hostDstS.setVal(ghostVal);
    hostSrc.copyTo(hostDstL);
    hostSrc.copyTo(hostDstS);
    EXPECT_TRUE(testCopyTo(hostSrc, hostDstL));
    EXPECT_TRUE(testCopyTo(hostSrc, hostDstS));
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
