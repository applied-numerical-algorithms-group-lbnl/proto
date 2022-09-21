#include <gtest/gtest.h>
#include "Proto.H"
#include "Lambdas.H"

using namespace Proto;
using namespace Operator;

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


TEST(Operator, Convolve) {
    int domainSize = 64;
    double offset = 0.1;
    int numIter = 3;
    double errTol = 1e-5;
    double rateTol = 0.1;

    HDF5Handler h5;
    double errNorm[numIter];
    for (int n = 0; n < numIter; n++)
    {
        double dx = 1.0/domainSize;
        Box B = Box::Cube(domainSize);
        BoxData<double, 1> data(B.grow(1));
        BoxData<double, 1> soln(B);
        BoxData<double, 1> err(B,0);
        forallInPlace_p(f_phi,     data, dx, offset);
        forallInPlace_p(f_phi_avg, soln, dx, offset);
        auto test = convolve(data);
        err += test;
        err -= soln;
        EXPECT_EQ(test.box(), B);
        errNorm[n] = err.absMax();
        EXPECT_LT(errNorm[n], errTol);
        domainSize *= 2; 
    }
    for (int ii = 1; ii < numIter; ii++)
    {
        double rate = log(errNorm[ii-1]/errNorm[ii])/log(2.0);
        double rateErr = abs(4-rate);
        EXPECT_LT(rateErr, rateTol);
    }
}

TEST(Operator, Deconvolve) {
    int domainSize = 64;
    double offset = 0.1;
    int numIter = 3;
    double errTol = 1e-5;
    double rateTol = 0.1;

    double errNorm[numIter];
    for (int n = 0; n < numIter; n++)
    {
        double dx = 1.0/domainSize;
        Box B = Box::Cube(domainSize);
        BoxData<double, 1> data(B.grow(1));
        BoxData<double, 1> soln(B);
        BoxData<double, 1> err(B,0);
        forallInPlace_p(f_phi_avg,  data, dx, offset);
        forallInPlace_p(f_phi,      soln, dx, offset);
        auto test = deconvolve(data);
        err += test;
        err -= soln;
        EXPECT_EQ(test.box(), B);
        errNorm[n] = err.absMax();
        EXPECT_LT(errNorm[n], errTol);
        domainSize *= 2; 
    }

    for (int ii = 1; ii < numIter; ii++)
    {
        double rate = log(errNorm[ii-1]/errNorm[ii])/log(2.0);
        double rateErr = abs(4-rate);
        EXPECT_LT(rateErr, rateTol);
    }
}

TEST(Operator, InitConvolve)
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
        Operator::initConvolve(hostData, f_phi, dx, offset);
        soln.initialize(f_phi_avg, dx, offset);
        hostErr[nn] = 0;
#ifdef PROTO_CUDA 
        LevelBoxData<double, 1, DEVICE> deviData(layout, Point::Ones());
        Operator::initConvolve(deviData, f_phi, dx, offset);
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

TEST(Operator, Cofactor)
{
    HDF5Handler h5;
    int domainSize = 32;
    double dx = 1.0/domainSize;
    Box domainBox = Box::Cube(domainSize+1);
    BoxData<double, DIM, HOST> X(domainBox);
    forallInPlace_p(f_iotaCorner, X, dx);

    std::array<BoxData<double, DIM, HOST>, DIM> N;
    for (int dir = 0; dir < DIM; dir++)
    {
        N[dir] = Operator::_cofactor(X, dir);
        h5.writePatch(dx, N[dir], "N%i", dir);
    }
    std::cout << "X: " << X.box() << " | N0: " << N[0].box() << " | N1: " << N[1].box() << std::endl;
    h5.writePatch(dx, X, "X");
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
