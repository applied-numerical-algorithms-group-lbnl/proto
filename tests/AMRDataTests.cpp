#include <gtest/gtest.h>
#include "Proto.H"
#include "Lambdas.H"

using namespace Proto;


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

TEST(AMRData, Initialize) {
    int domainSize = 32;
    int numLevels = 3;
    Point offset(1,2,3,4,5,6);
    Point k(1,2,3,4,5,6);
    double dx = 1.0/domainSize;
    Point refRatio = Point::Ones(2);
    Point boxSize = Point::Ones(16);
    auto grid = telescopingGrid(domainSize, numLevels, refRatio, boxSize);
    AMRData<double, 1, HOST> hostData(grid, Point::Ones());
    hostData.initialize(dx, f_phi, k, offset);
#ifdef PROTO_CUDA
    AMRData<double, 1, DEVICE> deviData(grid, Point::Ones());
    deviData.initialize(dx, f_phi, k, offset);
#endif
    for (int lvl = 0; lvl < numLevels; lvl++)
    {
        double dx_lvl = dx / pow(refRatio[0], lvl);
        for (auto iter : grid[lvl])
        {
            auto& hostData_i = hostData[lvl][iter];
            int N = hostData_i.size();
            Box B = hostData_i.box();
            BoxData<double, 1, HOST> soln_i(B);
            forallInPlace_p(f_phi, soln_i, dx_lvl, k, offset);
            EXPECT_TRUE(compareBoxData(soln_i, hostData_i));
#ifdef PROTO_CUDA
            BoxData<double, 1, HOST> tmpData_i(B);
            auto& deviData_i = deviData[lvl][iter];
            deviData_i.copyTo(tmpData_i);
            EXPECT_TRUE(compareBoxData(soln_i, tmpData_i));
#endif
        }
    }
}

TEST(AMRData, InitConvolve)
{
    int domainSize = 32;
    int numLevels = 3;
    Point offset(1,2,3,4,5,6);
    Point k(1,2,3,4,5,6);
    int numIter = 3;
    Point refRatio = Point::Ones(2);
    Point boxSize = Point::Ones(16);
    double hostErr[numIter];
#ifdef PROTO_CUDA
    double deviErr[numIter];
#endif
    for (int nn = 0; nn < numIter; nn++)
    {
        double dx = 1.0/domainSize;
        auto grid = telescopingGrid(domainSize, numLevels, refRatio, boxSize);
        AMRData<double, 1, HOST> hostData(grid, Point::Ones());
        AMRData<double, 1, HOST> soln(grid, Point::Ones());
        AMRData<double, 1, HOST> error(grid, Point::Ones());
        Operator::initConvolve(hostData, dx, f_phi, k, offset);
        soln.initialize(dx, f_phi_avg, k, offset);
        hostErr[nn] = 0;
#ifdef PROTO_CUDA 
        AMRData<double, 1, DEVICE> deviData(grid, Point::Ones());
        Operator::initConvolve(deviData, dx, f_phi, k, offset);
        deviErr[nn] = 0;
#endif
        for (int lvl = 0; lvl < numLevels; lvl++)
        {
            for (auto iter : grid[lvl])
            {
                auto& hostData_i = hostData[lvl][iter];
                auto& soln_i = soln[lvl][iter];
                auto& error_i = error[lvl][iter];
                hostData_i.copyTo(error_i);
                error_i -= soln_i;
            }
            hostErr[nn] = std::max(hostErr[nn], error[lvl].absMax());
#ifdef PROTO_CUDA
            for (auto iter : grid[lvl])
            {
                auto& deviData_i = deviData[lvl][iter];
                auto& soln_i = soln[lvl][iter];
                auto& error_i = error[lvl][iter];
                deviData_i.copyTo(error_i);
                error_i -= soln_i;
            }
            deviErr[nn] = std::max(deviErr[nn], error[lvl].absMax());
#endif
        }
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
