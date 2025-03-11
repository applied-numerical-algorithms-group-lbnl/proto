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

TEST(Operator, Convolve) {
    int domainSize = 64;
    Point offset(1,2,3,4,5,6);
    Point k(1,2,3,4,5,6);
    int numIter = 3;
    double errTol = 1e-5;
    double rateTol = 0.1;

    {
        ConvolveOp<double> convolveOp;
        int size2 = 8;
        int size4 = 10;
        int rangeSize = 8;
        Box domain2ndOrder = Box::Cube(size2);
        Box domain4thOrder = Box::Cube(size4);
        auto rangeBox = convolveOp.getRange(domain4thOrder, domain2ndOrder);
        EXPECT_EQ(rangeBox, domain2ndOrder.grow(-1));
        auto domainBoxes = convolveOp.domains(Box::Cube(rangeSize));
        EXPECT_EQ(domainBoxes.size(), 2);
        EXPECT_EQ(domainBoxes[0], Box::Cube(rangeSize));
        EXPECT_EQ(domainBoxes[1], Box::Cube(rangeSize).grow(1));
    }

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

    for (int dd = 0; dd < DIM; dd++)
    {
        int s2 = 8;
        int s4 = 10;
        ConvolveFaceOp<double> op(dd);
        std::vector<Box> domains;
        Box b2 = Box::Cube(s2);
        Box b4 = Box::Cube(s4);
        domains.push_back(b4);
        domains.push_back(b2);
        Box range = op.range(domains);
        Point ghost = Point::Ones() - Point::Basis(dd);
        EXPECT_EQ(range, b2.grow(-ghost));
        auto computedDomains = op.domains(range);
        EXPECT_EQ(computedDomains[0], range);
        EXPECT_EQ(computedDomains[1], range.grow(ghost));
    }

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
        
            auto side = Side::Hi;
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

    {
        DeconvolveOp<double> op;
        int size2 = 8;
        int size4 = 10;
        int rangeSize = 8;
        Box domain2ndOrder = Box::Cube(size2);
        Box domain4thOrder = Box::Cube(size4);
        auto rangeBox = op.getRange(domain4thOrder, domain2ndOrder);
        EXPECT_EQ(rangeBox, domain2ndOrder.grow(-1));
        auto domainBoxes = op.domains(Box::Cube(rangeSize));
        EXPECT_EQ(domainBoxes.size(), 2);
        EXPECT_EQ(domainBoxes[0], Box::Cube(rangeSize));
        EXPECT_EQ(domainBoxes[1], Box::Cube(rangeSize).grow(1));
    }

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

    for (int dd = 0; dd < DIM; dd++)
    {
        int s2 = 8;
        int s4 = 10;
        DeconvolveFaceOp<double> op(dd);
        std::vector<Box> domains;
        Box b2 = Box::Cube(s2);
        Box b4 = Box::Cube(s4);
        domains.push_back(b4);
        domains.push_back(b2);
        Box range = op.range(domains);
        Point ghost = Point::Ones() - Point::Basis(dd);
        EXPECT_EQ(range, b2.grow(-ghost));
        auto computedDomains = op.domains(range);
        EXPECT_EQ(computedDomains[0], range);
        EXPECT_EQ(computedDomains[1], range.grow(ghost));
    }

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

TEST(Operator, FaceDiff)
{
    #if PR_VERBOSE > 0
    HDF5Handler h5;
    #endif

    int domainSize = 32;
    Box rangeBox = Box::Cube(domainSize);
    Box domainBox;
    std::vector<FaceDiff<double>> faceDiff;
    for (int dd = 0; dd < DIM; dd++)
    {
        faceDiff.push_back(FaceDiff<double>(dd, 4));
        Box domain_i = faceDiff[dd].domains(rangeBox)[0];
        domainBox += domain_i;
    }
    

    int numIter = 3;
    double err[numIter];
    for (int nn = 0; nn < numIter; nn++)
    {
        err[nn] = 0.0;
        double dx = 1.0/domainSize;
        BoxData<double, 2> phi(domainBox);

        Point offset(1,2,3,4,5,6);
        Point k(1,2,3,4,5,6);
        forallInPlace_p(f_phi_avg,  phi, dx, k, offset);

        for (int norm = 0; norm < DIM; norm++)
        {
            BoxData<double, 2, HOST, DIM> gradPhi = faceDiff[norm](phi);
            for (int dd = 0; dd < DIM; dd++)
            {
                BoxData<double, 2> D0(rangeBox);
                forallInPlace_p(f_Dphi_avg_face,  D0, dx, k, offset, dd, norm);
                BoxData<double, 2> gradComponent = plane(gradPhi, dd);
                gradComponent /= dx;
                BoxData<double, 2> errorComponent(rangeBox);
                gradComponent.copyTo(errorComponent);
                errorComponent -= D0;
                #if PR_VERBOSE > 0
                h5.writePatch(D0, "SolnDPhiDX%i_N%i", dd, norm);
                h5.writePatch(gradComponent, "TestDPhiDX%i_N%i", dd, norm);
                h5.writePatch(errorComponent, "ErrDPhiDX%i_N%i", dd, norm);
                #endif
                err[nn] = max(err[nn], errorComponent.absMax());
            }
        }
        domainSize *= 2;
    }
    for (int ii = 1; ii < numIter; ii++)
    {
        double rate = log(err[ii-1]/err[ii])/log(2.0);
        #if PR_VERBOSE > 0
        std::cout << "Error (Max Norm): " << err[ii] << " | Convergence Rate: " << rate << std::endl;
        #endif
        EXPECT_NEAR(rate, 4.0, 0.1);
    }
}

TEST(Operator, FaceAverageProduct)
{
    #if PR_VERBOSE > 0
    HDF5Handler h5;
    #endif

    constexpr int C = 2;
    int domainSize = 32;
    int numIter = 3;
    double error[numIter];
    for (int nn = 0; nn < numIter; nn++)
    {
        double dx = 1.0/domainSize;
        Box rangeBox = Box::Cube(domainSize);
        Point offset_A(1,2,3,4,5,6);
        Point offset_B(2,1,4,3,6,5);
        Point k_A(1,2,3,4,5,6);
        Point k_B(2,1,4,3,6,5);
        auto side = Side::Lo;
        error[nn] = 0;
        for (int normDir = 0; normDir < DIM; normDir++)
        {
            auto productOp = FaceAverageProductOp<double>(normDir);
            auto domainBoxes = productOp.domains(rangeBox);
            
            auto AVG = ConvolveFaceOp<double>(normDir);

            Box avgDomain = AVG.domainUnion(rangeBox);
            BoxData<double, C> A(avgDomain);
            BoxData<double, C> B(avgDomain);
            BoxData<double, C> AB(avgDomain);
            BoxData<double, C> solnData(rangeBox);

            
            forallInPlace_p(f_phi_face,  A, dx, k_A, offset_A, normDir, side);
            forallInPlace_p(f_phi_face,  B, dx, k_B, offset_B, normDir, side);

            AB.setVal(0);
            AB.incrementProduct(A, B);
            AVG(solnData, AB, AB);

            BoxData<double, C> A4_avg(domainBoxes[0]);
            BoxData<double, C> B4_avg(domainBoxes[1]);
            BoxData<double, C> A2_avg(domainBoxes[2]);
            BoxData<double, C> B2_avg(domainBoxes[3]);
            BoxData<double, C> testData(rangeBox);

            forallInPlace_p(f_phi_face_avg,  A4_avg, dx, k_A, offset_A, normDir, side);
            forallInPlace_p(f_phi_face_avg,  B4_avg, dx, k_B, offset_B, normDir, side);
            forallInPlace_p(f_phi_face_avg,  A2_avg, dx, k_A, offset_A, normDir, side);
            forallInPlace_p(f_phi_face_avg,  B2_avg, dx, k_B, offset_B, normDir, side);

            productOp(testData, A4_avg, B4_avg, A2_avg, B2_avg);
            
            BoxData<double, C> errData(rangeBox);
            errData.setVal(0);
            errData += testData;
            errData -= solnData;

            #if PR_VERBOSE > 0
            h5.writePatch(testData, "TEST_PROD_FACE_DAT_D%i_N%i", normDir, nn);
            h5.writePatch(A, "TEST_PROD_FACE_A_D%i_N%i", normDir, nn);
            h5.writePatch(B, "TEST_PROD_FACE_B_D%i_N%i", normDir, nn);
            h5.writePatch(AB, "TEST_PROD_FACE_AB_D%i_N%i", normDir, nn);
            h5.writePatch(solnData, "TEST_PROD_FACE_SLN_D%i_N%i", normDir, nn);
            h5.writePatch(errData, "TEST_PROD_FACE_ERR_D%i_N%i", normDir, nn);
            #endif

            error[nn] = max(error[nn], errData.absMax());
        }
        domainSize *= 2;
    }
    for (int nn = 1; nn < numIter; nn++)
    {
        double rate = log(error[nn-1]/error[nn])/ log(2.0);
        #if PR_VERBOSE > 0
        std::cout << "Error (Max Norm): " << error[nn] << " | Convergence Rate: " << rate << std::endl;
        #endif
        EXPECT_NEAR(rate, 4.0, 0.25);
    }
}

#if 0
#if DIM==3
#ifdef PR_MMB
TEST(Operator, Cofactor)
{
#ifdef PR_HDF5
    HDF5Handler h5;
#endif
    int domainSize_0 = 16;
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
            FluxBoxData<double, DIM, HOST> NT(domainBox);
            Array<BoxData<double, DIM, HOST>, DIM> NT0;
            for (int dir = 0; dir < DIM; dir++)
            {
                NT[dir] = Operator::cofactor(X, dir);
            }
            BoxData<double, 1, HOST> J;
            J = Operator::jacobian(X, NT);
            //J /= dv;

            #if PR_VERBOSE > 0
            pr_out() << "Coordinate Values: " << std::endl;
            X.printData();
            pr_out() << "Cofactors: " << std::endl;
            for (int dd = 0; dd < DIM; dd++)
            {
                NT[dd].printData();
            }
            pr_out() << "Jacobian: " << std::endl;
            J.printData();
            #endif

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
