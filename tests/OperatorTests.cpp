#include <gtest/gtest.h>
#include "Proto.H"
#include "TestFunctions.H"

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
    Array<bool, DIM> periodicity;
    periodicity.fill(true);
    ProblemDomain domain(domainBox, periodicity);
    return DisjointBoxLayout(domain, patches, boxSize);
}
#if 0
TEST(Operator, Convolve) {
    int domainSize = 64;
    Point offset(1,2,3,4,5,6);
    Point k(1,2,3,4,5,6);
    int numIter = 3;
    double errTol = 1e-5;
    double rateTol = 0.1;

#ifdef PR_HDF5
    HDF5Handler h5;
#endif
    double errNorm[numIter];
    for (int n = 0; n < numIter; n++)
    {
        double dx = 1.0/domainSize;
        Box B = Box::Cube(domainSize);
        BoxData<double, 1> data(B.grow(1));
        BoxData<double, 1> soln(B);
        BoxData<double, 1> err(B,0);
        forallInPlace_p(f_phi,     data, dx, k, offset);
        forallInPlace_p(f_phi_avg, soln, dx, k, offset);
        auto test = convolve(data);
        err += test;
        err -= soln;
        EXPECT_EQ(test.box(), B);
        errNorm[n] = err.absMax();
        //EXPECT_LT(errNorm[n], errTol); TODO: figure out a proper error bound
        domainSize *= 2; 
    }
    for (int ii = 1; ii < numIter; ii++)
    {
        double rate = log(errNorm[ii-1]/errNorm[ii])/log(2.0);
        double rateErr = abs(4.0-rate);
        EXPECT_LT(rateErr, rateTol);
    }
}

TEST(Operator, ConvolveFace) {
    int domainSize = 64;
    Point offset(1,2,3,4,5,6);
    Point k(1,2,3,4,5,6);
    int numIter = 3;
    double errTol = 1e-5;
    double rateTol = 0.1;

#ifdef PR_HDF5
    HDF5Handler h5;
#endif
    double errNorm[numIter][DIM];
    for (int n = 0; n < numIter; n++)
    {
        double dx = 1.0/domainSize;
        Box B = Box::Cube(domainSize);
        Array<BoxData<double, 1, HOST>, DIM> data;
        Array<BoxData<double, 1, HOST>, DIM> soln;
        Array<BoxData<double, 1, HOST>, DIM> err;
        for (int dir = 0; dir < DIM; dir++)
        {
            Box Bd = B.grow(1).grow(dir,-1);
            data[dir].define(Bd);
            soln[dir].define(B);
            err[dir].define(B);
        
            forallInPlace_p(f_phi_face, data[dir], dx, k, offset, dir);
            forallInPlace_p(f_phi_face_avg, soln[dir], dx, k, offset, dir);
            err[dir].setVal(0);

            auto test = convolveFace(data[dir], dir);
            err[dir] += test;
            err[dir] -= soln[dir];
            EXPECT_EQ(test.box(), B);
            errNorm[n][dir] = err[dir].absMax();
        }
        domainSize *= 2; 
    }
    for (int ii = 1; ii < numIter; ii++)
    {
        for (int dir = 0; dir < DIM; dir++)
        {
            double rate = log(errNorm[ii-1][dir]/errNorm[ii][dir])/log(2.0);
            double rateErr = abs(4.0-rate);
            EXPECT_LT(rateErr, rateTol);
        }
    }
}


TEST(Operator, Deconvolve) {
    int domainSize = 64;
    Point offset(1,2,3,4,5,6);
    Point k(1,2,3,4,5,6);
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
        forallInPlace_p(f_phi_avg,  data, dx, k, offset);
        forallInPlace_p(f_phi,      soln, dx, k, offset);
        auto test = deconvolve(data);
        err += test;
        err -= soln;
        EXPECT_EQ(test.box(), B);
        errNorm[n] = err.absMax();
        //EXPECT_LT(errNorm[n], errTol); TODO: figure out a proper error bound
        domainSize *= 2; 
    }

    for (int ii = 1; ii < numIter; ii++)
    {
        double rate = log(errNorm[ii-1]/errNorm[ii])/log(2.0);
        double rateErr = abs(4-rate);
        EXPECT_LT(rateErr, rateTol);
    }
}

TEST(Operator, DeconvolveFace) {
    int domainSize = 64;
    Point offset(1,2,3,4,5,6);
    Point k(1,2,3,4,5,6);
    int numIter = 3;
    double errTol = 1e-5;
    double rateTol = 0.1;

#ifdef PR_HDF5
    HDF5Handler h5;
#endif
    double errNorm[numIter][DIM];
    for (int n = 0; n < numIter; n++)
    {
        double dx = 1.0/domainSize;
        Box B = Box::Cube(domainSize);
        Array<BoxData<double, 1, HOST>, DIM> data;
        Array<BoxData<double, 1, HOST>, DIM> soln;
        Array<BoxData<double, 1, HOST>, DIM> err;
        for (int dir = 0; dir < DIM; dir++)
        {
            Box Bd = B.grow(1).grow(dir,-1);
            data[dir].define(Bd);
            soln[dir].define(B);
            err[dir].define(B);
        
            forallInPlace_p(f_phi_face_avg, data[dir], dx, k, offset, dir);
            forallInPlace_p(f_phi_face, soln[dir], dx, k, offset, dir);
            err[dir].setVal(0);

            auto test = deconvolveFace(data[dir], dir);
            err[dir] += test;
            err[dir] -= soln[dir];
            EXPECT_EQ(test.box(), B);
            errNorm[n][dir] = err[dir].absMax();
        }
        domainSize *= 2; 
    }
    for (int ii = 1; ii < numIter; ii++)
    {
        for (int dir = 0; dir < DIM; dir++)
        {
            double rate = log(errNorm[ii-1][dir]/errNorm[ii][dir])/log(2.0);
            double rateErr = abs(4.0-rate);
            EXPECT_LT(rateErr, rateTol);
        }
    }
}


TEST(Operator, InitConvolve)
{
    int domainSize = 32;
    Point offset(1,2,3,4,5,6);
    Point k(1,2,3,4,5,6);
    int numIter = 3;
    Point boxSize = Point::Ones(16);
    double hostErr[numIter];
#ifdef PROTO_ACCEL
    double deviErr[numIter];
#endif
    for (int nn = 0; nn < numIter; nn++)
    {
        double dx = 1.0/domainSize;
        auto layout = testLayout(domainSize, boxSize);
        
        LevelBoxData<double, 1, HOST> hostData(layout, Point::Ones());
        LevelBoxData<double, 1, HOST> soln(layout, Point::Ones());
        LevelBoxData<double, 1, HOST> error(layout, Point::Ones());
        Operator::initConvolve(hostData, f_phi, dx, k, offset);
        soln.initialize(f_phi_avg, dx, k, offset);
        hostErr[nn] = 0;
#ifdef PROTO_ACCEL 
        LevelBoxData<double, 1, DEVICE> deviData(layout, Point::Ones());
        Operator::initConvolve(deviData, f_phi, dx, k, offset);
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
#ifdef PROTO_ACCEL
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
#ifdef PROTO_ACCEL
        double deviRate_i = log(deviErr[ii-1]/deviErr[ii])/log(2.0);
        EXPECT_TRUE(abs(rate - deviRate_i) < 0.1);
#endif
    }
}
#endif
#if DIM==3
#ifdef PR_MMB
TEST(Operator, Cofactor)
{
#ifdef PR_HDF5
    HDF5Handler h5;
#endif
    int domainSize_0 = 32;
    int thickness = 1;
    double r0 = 0.9;
    double r1 = 1.1;
    int ghostSize = 3;

    int N = 3;
    for (int tt = 0; tt < 2; tt++)
    {
        //tt == 0 -> classical sphere
        //tt == 1 -> cube sphere
        int domainSize = domainSize_0;
        double dv0; 
        switch (tt)
        {
            case 0: dv0 = 2.0*M_PI*M_PI/2.0*(r1-r0); break;
            case 1: dv0 = M_PI*M_PI/4.0*(r1-r0); break;
        }
        //dv0 = 1.0;
        double JErrNorm[N];
        for (int nn = 0; nn < N; nn++)
        {
            Array<double, DIM> dx = 1.0/domainSize;
            dx[0] = 1.0/thickness;
            double dv = dv0;
            for (int dir = 0; dir < DIM; dir++) { dv *= dx[dir]; }

            // Define Mapping
            Point domainBoxSize(thickness, domainSize, domainSize);
            Box domainBox(domainBoxSize);
            domainBox = domainBox.grow(ghostSize+1);
            BoxData<double, DIM, HOST> X(domainBox);
            BoxData<double, 1, HOST> J0(domainBox);
            switch(tt)
            {
                case 0: forallInPlace_p(f_classicSphereMap, X, J0, dx, r0, r1); break;
                //case 1: forallInPlace_p(f_cubeSphereMap, X, J0, dx, r0, r1); break;
            }

            // Compute Metrics
            FluxBoxData<double, DIM, HOST> NT;
            Array<BoxData<double, DIM, HOST>, DIM> NT0;
            for (int dir = 0; dir < DIM; dir++)
            {
                NT[dir] = Operator::cofactor(X, dir);
            }
            BoxData<double, 1, HOST> J;
            J = Operator::jacobian(X, NT);
            //J /= dv;

            auto JAvg = Operator::convolve(J0);
            JAvg *= dv;
            BoxData<double, 1> JAvg2(J.box());
            JAvg.copyTo(JAvg2);
            BoxData<double, 1> JErr(J.box());
            J.copyTo(JErr);
            JErr -= JAvg;
            JErrNorm[nn] = JErr.sumAbs()*dv;
#if PR_VERBOSE > 0
            std::cout << "Max Norm Error in J: " << JErrNorm[nn] << std::endl;
            h5.writePatch({"x", "y", "z"}, dx, X, "X_T%i_N%i", tt, nn);
            h5.writePatch({"J"}, dx, JAvg2, "J0_T%i_N%i", tt, nn);
            h5.writePatch({"J"}, dx, J, "J_T%i_N%i", tt, nn);
            h5.writePatch({"Err"}, dx, JErr, "JErr_T%i_N%i", tt, nn);
#endif
            domainSize *= 2;
        }
        
        for (int ii = 1; ii < N; ii++)
        {
            double rate = log(JErrNorm[ii-1]/JErrNorm[ii])/log(2.0);
            double rateErr = std::abs(rate - 4);
#if PR_VERBOSE > 0
            std::cout << "Convergence Rate: " << rate << std::endl;
#endif
            EXPECT_LT(rateErr, 0.3);
        }

    }
}
#endif
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
