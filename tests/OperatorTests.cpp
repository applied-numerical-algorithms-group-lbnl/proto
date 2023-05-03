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

#ifdef DIM==3
PROTO_KERNEL_START
void f_cubeSphereMapF(Point& a_pt, Var<double,3>& a_X, Var<double,1>& a_J, Array<double, 3> a_dx, double a_r0, double a_r1)
{
    // subscript 0 are at corners
    Array<double, DIM> x  = a_pt;
    x += 0.5;
    x *= a_dx;
    Array<double, DIM> x0 = x - (0.5*a_dx);

    double r0   = a_r0 + (a_r1-a_r0)*x0[0];
    double xi0  = -M_PI/4.0 + M_PI/2.0*x0[1];
    double eta0 = -M_PI/4.0 + M_PI/2.0*x0[2];
    
    double r   = a_r0 + (a_r1-a_r0)*x[0];
    double xi  = -M_PI/4.0 + M_PI/2.0*x[1];
    double eta = -M_PI/4.0 + M_PI/2.0*x[2];

    double X0 = tan(xi0);
    double Y0 = tan(eta0);
    double d0 = sqrt(1+X0*X0+Y0*Y0);
    
    double X = tan(xi);
    double Y = tan(eta);
    double d = sqrt(1+X*X+Y*Y);

    a_J(0) = r*r*(1+X*X)*(1+Y*Y)/(d*d*d);
    a_X(0) = r0/d0;
    a_X(1) = r0*X0/d0;
    a_X(2) = r0*Y0/d0;
}
PROTO_KERNEL_END(f_cubeSphereMapF, f_cubeSphereMap);

PROTO_KERNEL_START
void f_classicSphereMapF(Point& a_pt, Var<double,3>& a_X, Var<double,1>& a_J, Array<double, 3> a_dx, double a_r0, double a_r1)
{
    Array<double, DIM> x = a_pt;
    x += 0.5;
    x *= a_dx;
    Array<double, DIM> x0 = x - (0.5*a_dx);
    Array<double, DIM> x1 = x0 + a_dx;

    double r      = a_r0 + (a_r1-a_r0)*x[0];
    double phi    = M_PI/4.0 + M_PI/2.0*x[1];
    double theta  = 2.0*M_PI*x[2];
    
    double r0      = a_r0 + (a_r1-a_r0)*x0[0];
    double phi0    = M_PI/4.0 + M_PI/2.0*x0[1];
    double theta0  = 2.0*M_PI*x0[2];
    
    double r1      = a_r0 + (a_r1-a_r0)*x1[0];
    double phi1    = M_PI/4.0 + M_PI/2.0*x1[1];
    double theta1  = 2.0*M_PI*x1[2];
    
    a_J(0) = r*r*sin(phi);
    
    a_X(0) = r0*cos(theta0)*sin(phi0);
    a_X(1) = r0*sin(theta0)*sin(phi0);
    a_X(2) = r0*cos(phi0);

    //a_N(0,0) =  r0*sin(phi)*sin(theta);
    //a_N(1,0) = -r0*sin(phi)*cos(theta);
    //a_N(2,0) =  r0*cos(phi)*cos(theta);
    //a_N(0,1) = -r*cos(phi0)*cos(theta);
    //a_N(1,1) =  r*cos(phi0)*sin(theta);
    //a_N(2,1) =  0;
    //a_N(0,2) =  -r*sin(phi)*cos(theta0);
    //a_N(1,2) =  -r*sin(phi)*sin(theta0);
    //a_N(2,2) =  -r*sin(phi);
}
PROTO_KERNEL_END(f_classicSphereMapF, f_classicSphereMap);

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
                case 1: forallInPlace_p(f_cubeSphereMap, X, J0, dx, r0, r1); break;
            }

            // Compute Metrics
            Array<BoxData<double, DIM, HOST>, DIM> NT;
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
#if PR_VERBOSE > 0
            std::cout << "Convergence Rate: " << rate << std::endl;
#endif
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
