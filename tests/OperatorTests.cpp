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
        #if PR_VERBOSE > 0
        std::cout << "Error (Max Norm): " << errNorm[ii] << " | Convergence Rate: " << rate << std::endl;
        #endif
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
        std::array<Box,2> domains;
        Box b2 = Box::Cube(s2);
        Box b4 = Box::Cube(s4);
        domains[0]=b4;
        domains[1]=b2;
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
            #if PR_VERBOSE > 0
            std::cout << "Dir: " << dir << " | Error (Max Norm): " << errNorm[ii][dir] << " | Convergence Rate: " << rate << std::endl;
            #endif
            EXPECT_LT(rateErr, rateTol);
        }
    }
}

TEST(Operator, ConvolveEdge) {
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
        ConvolveEdgeOp<double> op(dd);
        std::array<Box,2> domains;
        Box b2 = Box::Cube(s2);
        Box b4 = Box::Cube(s4);
        domains[0]=b4;
        domains[1]=b2;
        Box range = op.range(domains);
        Point ghost = Point::Basis(dd);
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
            ConvolveEdgeOp<double> convolveEdge(dir);
            Box Bd = B.grow(dir,1);
            data[dir].define(Bd);
            soln[dir].define(B);
            err[dir].define(B);
        
            forallInPlace_p(f_phi_edge, data[dir], dx, k, offset, dir);
            forallInPlace_p(f_phi_edge_avg, soln[dir], dx, k, offset, dir);
            err[dir].setVal(0);

            auto test = convolveEdge(data[dir], data[dir]);
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
            #if PR_VERBOSE > 0
            std::cout << "Dir: " << dir << " | Error (Max Norm): " << errNorm[ii][dir] << " | Convergence Rate: " << rate << std::endl;
            #endif
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
        #if PR_VERBOSE > 0
        std::cout << "Error (Max Norm): " << errNorm[ii] << " | Convergence Rate: " << rate << std::endl;
        #endif
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
        std::array<Box,2> domains;
        Box b2 = Box::Cube(s2);
        Box b4 = Box::Cube(s4);
        domains[1] = b2;
        domains[0] = b4;
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
            #if PR_VERBOSE > 0
            std::cout << "Dir: " << dir << " | Error (Max Norm): " << errNorm[ii][dir] << " | Convergence Rate: " << rate << std::endl;
            #endif
            EXPECT_LT(rateErr, rateTol);
        }
    }
}

TEST(Operator, DeconvolveEdge) {
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
        DeconvolveEdgeOp<double> op(dd);
        std::array<Box,2> domains;
        Box b2 = Box::Cube(s2);
        Box b4 = Box::Cube(s4);
        domains[0]=b4;
        domains[1]=b2;
        Box range = op.range(domains);
        Point ghost = Point::Basis(dd);
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
            Box Bd = B.grow(dir,1);
            data[dir].define(Bd);
            soln[dir].define(B);
            err[dir].define(B);
        
            forallInPlace_p(f_phi_edge_avg, data[dir], dx, k, offset, dir);
            forallInPlace_p(f_phi_edge, soln[dir], dx, k, offset, dir);
            err[dir].setVal(0);

            DeconvolveEdgeOp<double> deconvolveEdge(dir);
            auto test = deconvolveEdge(data[dir], data[dir]);
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
            #if PR_VERBOSE > 0
            std::cout << "Dir: " << dir << " | Error (Max Norm): " << errNorm[ii][dir] << " | Convergence Rate: " << rate << std::endl;
            #endif
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

TEST(Operator, FaceAverageDiffOp)
{
    #if PR_VERBOSE > 0
    HDF5Handler h5;
    #endif

    int domainSize = 32;
    Box rangeBox = Box::Cube(domainSize);
    Box domainBox;
    std::vector<FaceAverageDiffOp<double>> faceDiff;
    for (int dd = 0; dd < DIM; dd++)
    {
        faceDiff.push_back(FaceAverageDiffOp<double>(dd, 4));
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

TEST(Operator, CellAverageProduct)
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

        CellAverageProductOp<double> productOp;
        auto domainBoxes = productOp.domains(rangeBox);
        
        ConvolveOp<double> AVG;

        Box avgDomain = AVG.domainUnion(rangeBox);
        BoxData<double, C> A(avgDomain);
        BoxData<double, C> B(avgDomain);
        BoxData<double, C> AB(avgDomain);
        BoxData<double, C> solnData(rangeBox);

        
        forallInPlace_p(f_phi,  A, dx, k_A, offset_A);
        forallInPlace_p(f_phi,  B, dx, k_B, offset_B);

        AB.setVal(0);
        AB.incrementProduct(A, B);
        AVG(solnData, AB, AB);

        BoxData<double, C> A4_avg(domainBoxes[0]);
        BoxData<double, C> B4_avg(domainBoxes[1]);
        BoxData<double, C> A2_avg(domainBoxes[2]);
        BoxData<double, C> B2_avg(domainBoxes[3]);
        BoxData<double, C> testData(rangeBox);

        forallInPlace_p(f_phi_avg,  A4_avg, dx, k_A, offset_A);
        forallInPlace_p(f_phi_avg,  B4_avg, dx, k_B, offset_B);
        forallInPlace_p(f_phi_avg,  A2_avg, dx, k_A, offset_A);
        forallInPlace_p(f_phi_avg,  B2_avg, dx, k_B, offset_B);

        productOp(testData, A4_avg, B4_avg, A2_avg, B2_avg);
        
        BoxData<double, C> errData(rangeBox);
        errData.setVal(0);
        errData += testData;
        errData -= solnData;

        #if PR_VERBOSE > 0
        h5.writePatch(testData, "TEST_PROD_CELL_DAT_N%i", nn);
        h5.writePatch(A, "TEST_PROD_CELL_A_N%i", nn);
        h5.writePatch(B, "TEST_PROD_CELL_B_N%i", nn);
        h5.writePatch(AB, "TEST_PROD_CELL_AB_N%i", nn);
        h5.writePatch(solnData, "TEST_PROD_CELL_SLN_N%i", nn);
        h5.writePatch(errData, "TEST_PROD_CELL_ERR_N%i", nn);
        #endif

        error[nn] = errData.absMax();
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

TEST(Operator, EdgeAverageProduct)
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
        for (int edgeDir = 0; edgeDir < DIM; edgeDir++)
        {
            auto productOp = EdgeAverageProductOp<double>(edgeDir);
            auto domainBoxes = productOp.domains(rangeBox);
            
            auto AVG = ConvolveEdgeOp<double>(edgeDir);

            Box avgDomain = AVG.domainUnion(rangeBox);
            BoxData<double, C> A(avgDomain);
            BoxData<double, C> B(avgDomain);
            BoxData<double, C> AB(avgDomain);
            BoxData<double, C> solnData(rangeBox);

            
            forallInPlace_p(f_phi_edge,  A, dx, k_A, offset_A, edgeDir);
            forallInPlace_p(f_phi_edge,  B, dx, k_B, offset_B, edgeDir);

            AB.setVal(0);
            AB.incrementProduct(A, B);
            AVG(solnData, AB, AB);

            BoxData<double, C> A4_avg(domainBoxes[0]);
            BoxData<double, C> B4_avg(domainBoxes[1]);
            BoxData<double, C> A2_avg(domainBoxes[2]);
            BoxData<double, C> B2_avg(domainBoxes[3]);
            BoxData<double, C> testData(rangeBox);

            forallInPlace_p(f_phi_edge_avg,  A4_avg, dx, k_A, offset_A, edgeDir);
            forallInPlace_p(f_phi_edge_avg,  B4_avg, dx, k_B, offset_B, edgeDir);
            forallInPlace_p(f_phi_edge_avg,  A2_avg, dx, k_A, offset_A, edgeDir);
            forallInPlace_p(f_phi_edge_avg,  B2_avg, dx, k_B, offset_B, edgeDir);

            productOp(testData, A4_avg, B4_avg, A2_avg, B2_avg);
            
            BoxData<double, C> errData(rangeBox);
            errData.setVal(0);
            errData += testData;
            errData -= solnData;

            #if PR_VERBOSE > 0
            h5.writePatch(testData, "TEST_PROD_EDGE_DAT_D%i_N%i", edgeDir, nn);
            h5.writePatch(A, "TEST_PROD_EDGE_A_D%i_N%i", edgeDir, nn);
            h5.writePatch(B, "TEST_PROD_EDGE_B_D%i_N%i", edgeDir, nn);
            h5.writePatch(AB, "TEST_PROD_EDGE_AB_D%i_N%i", edgeDir, nn);
            h5.writePatch(solnData, "TEST_PROD_EDGE_SLN_D%i_N%i", edgeDir, nn);
            h5.writePatch(errData, "TEST_PROD_EDGE_ERR_D%i_N%i", edgeDir, nn);
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

TEST(Operator, CellAverageQuotient)
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

        CellAverageQuotientOp<double> quotientOp;
        auto domainBoxes = quotientOp.domains(rangeBox);
        
        ConvolveOp<double> AVG;

        Box avgDomain = AVG.domainUnion(rangeBox);
        BoxData<double, C> A(avgDomain);
        BoxData<double, C> B(avgDomain);
        BoxData<double, C> AB(avgDomain);
        BoxData<double, C> solnData(rangeBox);
        
        forallInPlace_p(f_phi,  A, dx, k_A, offset_A);
        forallInPlace_p(f_phi,  B, dx, k_B, offset_B);
        B += 10.0;

        AB.setVal(0);
        A.copyTo(AB);
        AB /= B;
        AVG(solnData, AB, AB);

        BoxData<double, C> A4_avg(domainBoxes[0]);
        BoxData<double, C> B4_avg(domainBoxes[1]);
        BoxData<double, C> A2_avg(domainBoxes[2]);
        BoxData<double, C> B2_avg(domainBoxes[3]);
        BoxData<double, C> testData(rangeBox);

        forallInPlace_p(f_phi_avg,  A4_avg, dx, k_A, offset_A);
        forallInPlace_p(f_phi_avg,  B4_avg, dx, k_B, offset_B);
        forallInPlace_p(f_phi_avg,  A2_avg, dx, k_A, offset_A);
        forallInPlace_p(f_phi_avg,  B2_avg, dx, k_B, offset_B);
        B4_avg += 10.0;
        B2_avg += 10.0;

        quotientOp(testData, A4_avg, B4_avg, A2_avg, B2_avg);
        
        BoxData<double, C> errData(rangeBox);
        errData.setVal(0);
        errData += testData;
        errData -= solnData;

        #if PR_VERBOSE > 0
        h5.writePatch(testData, "TEST_QUOT_CELL_DAT_N%i", nn);
        h5.writePatch(A, "TEST_QUOT_CELL_A_N%i", nn);
        h5.writePatch(B, "TEST_QUOT_CELL_B_N%i", nn);
        h5.writePatch(AB, "TEST_QUOT_CELL_AB_N%i", nn);
        h5.writePatch(solnData, "TEST_QUOT_CELL_SLN_N%i", nn);
        h5.writePatch(errData, "TEST_QUOT_CELL_ERR_N%i", nn);
        #endif

        error[nn] = errData.absMax();
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

TEST(Operator, FaceAverageQuotient)
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
            auto quotientOp = FaceAverageQuotientOp<double>(normDir);
            auto domainBoxes = quotientOp.domains(rangeBox);
            
            auto AVG = ConvolveFaceOp<double>(normDir);

            Box avgDomain = AVG.domainUnion(rangeBox);
            BoxData<double, C> A(avgDomain);
            BoxData<double, C> B(avgDomain);
            BoxData<double, C> AB(avgDomain);
            BoxData<double, C> solnData(rangeBox);

            
            forallInPlace_p(f_phi_face,  A, dx, k_A, offset_A, normDir, side);
            forallInPlace_p(f_phi_face,  B, dx, k_B, offset_B, normDir, side);
            B += 10.0;

            AB.setVal(0);
            A.copyTo(AB);
            AB /= B;
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

            B4_avg += 10.0;
            B2_avg += 10.0;

            quotientOp(testData, A4_avg, B4_avg, A2_avg, B2_avg);
            
            BoxData<double, C> errData(rangeBox);
            errData.setVal(0);
            errData += testData;
            errData -= solnData;

            #if PR_VERBOSE > 0
            h5.writePatch(testData, "TEST_QUOT_FACE_DAT_D%i_N%i", normDir, nn);
            h5.writePatch(A, "TEST_QUOT_FACE_A_D%i_N%i", normDir, nn);
            h5.writePatch(B, "TEST_QUOT_FACE_B_D%i_N%i", normDir, nn);
            h5.writePatch(AB, "TEST_QUOT_FACE_AB_D%i_N%i", normDir, nn);
            h5.writePatch(solnData, "TEST_QUOT_FACE_SLN_D%i_N%i", normDir, nn);
            h5.writePatch(errData, "TEST_QUOT_FACE_ERR_D%i_N%i", normDir, nn);
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

TEST(Operator, CellAverageMatrixProduct)
{
    #if PR_VERBOSE > 0
    HDF5Handler h5;
    #endif

    constexpr int CL = 2;
    constexpr int DL = 4;
    constexpr int CR = 4;
    constexpr int DR = 1;
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

        CellAverageMatrixProductOp<double> matProductOp;
        auto domainBoxes = matProductOp.domains(rangeBox);
    
        BoxData<double, CL, HOST, DL> A4_avg(domainBoxes[0]);
        BoxData<double, CR, HOST, DR> B4_avg(domainBoxes[1]);
        BoxData<double, CL, HOST, DL> A2_avg(domainBoxes[2]);
        BoxData<double, CR, HOST, DR> B2_avg(domainBoxes[3]);

        forallInPlace_p(f_phi_avg,  A4_avg, dx, k_A, offset_A);
        forallInPlace_p(f_phi_avg,  B4_avg, dx, k_B, offset_B);
        forallInPlace_p(f_phi_avg,  A2_avg, dx, k_A, offset_A);
        forallInPlace_p(f_phi_avg,  B2_avg, dx, k_B, offset_B);

        auto testAB = matProductOp(A4_avg, B4_avg, A2_avg, B2_avg);
        auto solnAB = Operator::_cellMatrixProductAB(A4_avg, B4_avg, A2_avg, B2_avg);
        
        BoxData<double, CL, HOST, DR> errData(rangeBox);
        errData.setVal(0);
        errData += testAB;
        errData -= solnAB;

        #if PR_VERBOSE > 0
        h5.writePatch(testAB, "TEST_MAT_PROD_CELL_DAT_N%i", nn);
        h5.writePatch(solnAB, "TEST_MAT_PROD_CELL_SLN_N%i", nn);
        h5.writePatch(errData, "TEST_MAT_PROD_CELL_ERR_N%i", nn);
        #endif

        error[nn] = max(error[nn], errData.absMax());

        domainSize *= 2;
    }
    for (int nn = 1; nn < numIter; nn++)
    {
        double rate = log(error[nn-1]/error[nn])/ log(2.0);
        #if PR_VERBOSE > 0
        std::cout << "Error (w.r.t original implementation): " << error[nn] << " | Convergence Rate: " << rate << std::endl;
        #endif
        if (error[nn] > 1e-12)
        {
            EXPECT_NEAR(rate, 4.0, 0.25);
        }
    }
}

TEST(Operator, FaceAverageMatrixProduct)
{
    #if PR_VERBOSE > 0
    HDF5Handler h5;
    #endif

    constexpr int CL = 2;
    constexpr int DL = 4;
    constexpr int CR = 4;
    constexpr int DR = 1;
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
            auto matProductOp = FaceAverageMatrixProductOp<double>(normDir);
            auto domainBoxes = matProductOp.domains(rangeBox);
        
            BoxData<double, CL, HOST, DL> A4_avg(domainBoxes[0]);
            BoxData<double, CR, HOST, DR> B4_avg(domainBoxes[1]);
            BoxData<double, CL, HOST, DL> A2_avg(domainBoxes[2]);
            BoxData<double, CR, HOST, DR> B2_avg(domainBoxes[3]);

            BoxData<double, DL, HOST, CL> AT4_avg(domainBoxes[0]);
            BoxData<double, DR, HOST, CR> BT4_avg(domainBoxes[1]);
            BoxData<double, DL, HOST, CL> AT2_avg(domainBoxes[2]);
            BoxData<double, DR, HOST, CR> BT2_avg(domainBoxes[3]);

            forallInPlace_p(f_phi_face_avg,  A4_avg, dx, k_A, offset_A, normDir, side);
            forallInPlace_p(f_phi_face_avg,  B4_avg, dx, k_B, offset_B, normDir, side);
            forallInPlace_p(f_phi_face_avg,  A2_avg, dx, k_A, offset_A, normDir, side);
            forallInPlace_p(f_phi_face_avg,  B2_avg, dx, k_B, offset_B, normDir, side);
            
            forallInPlace_p(f_phi_face_avg,  AT4_avg, dx, k_A, offset_A, normDir, side);
            forallInPlace_p(f_phi_face_avg,  BT4_avg, dx, k_B, offset_B, normDir, side);
            forallInPlace_p(f_phi_face_avg,  AT2_avg, dx, k_A, offset_A, normDir, side);
            forallInPlace_p(f_phi_face_avg,  BT2_avg, dx, k_B, offset_B, normDir, side);

            auto testAB = matProductOp(A4_avg, B4_avg, A2_avg, B2_avg);
            auto solnAB = Operator::_faceMatrixProductAB(A4_avg, B4_avg, A2_avg, B2_avg, normDir);
            auto testATB = matProductOp.transposeLeft(AT4_avg, B4_avg, AT2_avg, B2_avg);
            auto solnATB = Operator::_faceMatrixProductATB(AT4_avg, B4_avg, AT2_avg, B2_avg, normDir);
            auto testABT = matProductOp.transposeRight(A4_avg, BT4_avg, A2_avg, BT2_avg);
            auto solnABT = Operator::_faceMatrixProductABT(A4_avg, BT4_avg, A2_avg, BT2_avg, normDir);
            
            BoxData<double, CL, HOST, DR> errAB(rangeBox);
            errAB.setVal(0);
            errAB += testAB;
            errAB -= solnAB;
            BoxData<double, CL, HOST, DR> errATB(rangeBox);
            errATB.setVal(0);
            errATB += testATB;
            errATB -= solnATB;
            BoxData<double, CL, HOST, DR> errABT(rangeBox);
            errABT.setVal(0);
            errABT += testABT;
            errABT -= solnABT;

            #if PR_VERBOSE > 0
            h5.writePatch(testAB, "TEST_MAT_PROD_FACE_AB_DAT_D%i_N%i", normDir, nn);
            h5.writePatch(solnAB, "TEST_MAT_PROD_FACE_AB_SLN_D%i_N%i", normDir, nn);
            h5.writePatch(errAB, "TEST_MAT_PROD_FACE_AB_ERR_D%i_N%i", normDir, nn);
            h5.writePatch(testATB, "TEST_MAT_PROD_FACE_ATB_DAT_D%i_N%i", normDir, nn);
            h5.writePatch(solnATB, "TEST_MAT_PROD_FACE_ATB_SLN_D%i_N%i", normDir, nn);
            h5.writePatch(errATB, "TEST_MAT_PROD_FACE_ATB_ERR_D%i_N%i", normDir, nn);
            h5.writePatch(testABT, "TEST_MAT_PROD_FACE_ABT_DAT_D%i_N%i", normDir, nn);
            h5.writePatch(solnABT, "TEST_MAT_PROD_FACE_ABT_SLN_D%i_N%i", normDir, nn);
            h5.writePatch(errABT, "TEST_MAT_PROD_FACE_ABT_ERR_D%i_N%i", normDir, nn);
            #endif

            error[nn] = max(error[nn], errAB.absMax());
            error[nn] = max(error[nn], errATB.absMax());
            error[nn] = max(error[nn], errABT.absMax());
        }
        domainSize *= 2;
    }
    for (int nn = 1; nn < numIter; nn++)
    {
        double rate = log(error[nn-1]/error[nn])/ log(2.0);
        #if PR_VERBOSE > 0
        std::cout << "Error (w.r.t original implementation): " << error[nn] << " | Convergence Rate: " << rate << std::endl;
        #endif
        if (error[nn] > 1e-12)
        {
            EXPECT_NEAR(rate, 4.0, 0.25);
        }
    }
}

TEST(Operator, FaceAverageTensorQuotient)
{
    #if PR_VERBOSE > 0
    HDF5Handler h5;
    #endif

    constexpr int C = 2;
    constexpr int D = 1;
    constexpr int E = 3;

    int domainSize = 32;

    double dx = 1.0/domainSize;
    Box rangeBox = Box::Cube(domainSize);
    Point offset_A(1,2,3,4,5,6);
    Point offset_B(2,1,4,3,6,5);
    Point k_A(1,2,3,4,5,6);
    Point k_B(2,1,4,3,6,5);
    auto side = Side::Lo;

    for (int normDir = 0; normDir < DIM; normDir++)
    {
        auto quotientOp = FaceAverageTensorQuotientOp<double>(normDir);
        auto domainBoxes = quotientOp.domains(rangeBox);

        BoxData<double, C, HOST, D, E> A4_avg(domainBoxes[0]);
        BoxData<double, 1> B4_avg(domainBoxes[1]);
        BoxData<double, C, HOST, D, E> A2_avg(domainBoxes[2]);
        BoxData<double, 1> B2_avg(domainBoxes[3]);
        BoxData<double, C, HOST, D, E> testData(rangeBox);

        forallInPlace_p(f_phi_face_avg,  A4_avg, dx, k_A, offset_A, normDir, side);
        forallInPlace_p(f_phi_face_avg,  B4_avg, dx, k_B, offset_B, normDir, side);
        forallInPlace_p(f_phi_face_avg,  A2_avg, dx, k_A, offset_A, normDir, side);
        forallInPlace_p(f_phi_face_avg,  B2_avg, dx, k_B, offset_B, normDir, side);

        B4_avg += 10.0;
        B2_avg += 10.0;

        quotientOp(testData, A4_avg, B4_avg, A2_avg, B2_avg);
        auto solnData = Operator::_faceTensorQuotient(A4_avg, B4_avg, A2_avg, B2_avg, normDir);
        solnData -= testData;
        EXPECT_LT(solnData.absMax(), 1e-12);
        
        auto QOp = FaceAverageQuotientOp<double>(normDir);

        for (int ee = 0; ee < E; ee++)
        for (int dd = 0; dd < D; dd++)
        for (int cc = 0; cc < C; cc++)
        {
            BoxData<double, 1> A4 = slice(A4_avg, cc, dd, ee);
            BoxData<double, 1> A2 = slice(A2_avg, cc, dd, ee);
            BoxData<double, 1> Q = slice(testData, cc, dd, ee);
            auto errData = QOp(A4, B4_avg, A2, B2_avg);
            errData -= Q;
            EXPECT_LT(errData.absMax(), 1e-12);
        }
    }
}

TEST(Operator, EdgeAverageCrossProduct)
{
    #if PR_VERBOSE > 0
    HDF5Handler h5;
    #endif

    int domainSize = 32;
    int numIter = 1;
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
        for (int edgeDir = 0; edgeDir < DIM; edgeDir++)
        {
            auto productOp = EdgeAverageCrossProductOp<double>(edgeDir);
            auto domainBoxes = productOp.domains(rangeBox);
            EXPECT_EQ(productOp.range(domainBoxes), rangeBox);

            BoxData<double, 3> A4_avg(domainBoxes[0]);
            BoxData<double, 3> B4_avg(domainBoxes[1]);
            BoxData<double, 3> A2_avg(domainBoxes[2]);
            BoxData<double, 3> B2_avg(domainBoxes[3]);
            BoxData<double, 3> testData(rangeBox);

            forallInPlace_p(f_phi_edge_avg,  A4_avg, dx, k_A, offset_A, edgeDir);
            forallInPlace_p(f_phi_edge_avg,  B4_avg, dx, k_B, offset_B, edgeDir);
            forallInPlace_p(f_phi_edge_avg,  A2_avg, dx, k_A, offset_A, edgeDir);
            forallInPlace_p(f_phi_edge_avg,  B2_avg, dx, k_B, offset_B, edgeDir);

            productOp(testData, A4_avg, B4_avg, A2_avg, B2_avg);
            
            auto solnData = Operator::_edgeCrossProduct3D(A4_avg, B4_avg, A2_avg, B2_avg, edgeDir);

            BoxData<double, 3> errData(rangeBox);
            errData.setVal(0);
            errData += testData;
            errData -= solnData;

            #if PR_VERBOSE > 0
            h5.writePatch(testData, "TEST_CROSS_PROD_EDGE_DAT_D%i_N%i", edgeDir, nn);
            h5.writePatch(solnData, "TEST_CROSS_PROD_EDGE_SLN_D%i_N%i", edgeDir, nn);
            h5.writePatch(errData, "TEST_CROSS_PROD_EDGE_ERR_D%i_N%i", edgeDir, nn);
            #endif

            error[nn] = max(error[nn], errData.absMax());
        }
        domainSize *= 2;
    }
    for (int nn = 1; nn < numIter; nn++)
    {
        double rate = log(error[nn-1]/error[nn])/ log(2.0);
        #if PR_VERBOSE > 0
        std::cout << "Error (w.r.t original implementation): " << error[nn] << " | Convergence Rate: " << rate << std::endl;
        #endif
        if (error[nn] > 1e-12)
        {
            EXPECT_NEAR(rate, 4.0, 0.25);
        }
    }
}
#endif
#if 1
namespace {
    
    PROTO_KERNEL_START
    void f_shearF(Var<double, DIM>& x, Var<double, DIM>& X, double stretch, double angle)
    {
        x(0) = stretch*X(0) + X(1)/tan(angle);
        for (int ii = 1; ii < DIM; ii++)
        {
            x(ii) = X(ii);
        }
    }
    PROTO_KERNEL_END(f_shearF, f_shear);
}

TEST(Operator, FaceAveragedSurfaceTransformOp)
{
    int domainSize = 32;
    double stretch = 2.0;
    double angle = M_PI/4.0;
    double dX = domainSize/domainSize;
    double dA = pow(dX,DIM);
    double a = stretch;
    double b = 1.0/tan(angle);
    Array<double, DIM> dx;
    dx.fill(dX);
    dx[0] *= stretch;
    dx[1] /= sin(angle);
    Array<double, DIM> da;

#if DIM == 2
    da[0] = dx[1];
    da[1] = dx[0];
#elif DIM == 3
    da[0] = dx[1]*dx[2];
    da[1] = dx[2]*dx[0];
    da[2] = dx[0]*dx[1]*sin(angle);
#endif

    FaceAveragedSurfaceTransformOp<double, DIM> adj0(0);
    FaceAveragedSurfaceTransformOp<double, DIM> adj1(1);

    Box range = Box::Cube(domainSize);

    auto [domain0] = adj0.domains(range);
    auto range0 = adj0.range({domain0});
    auto [domain1] = adj1.domains(range);
    auto range1 = adj1.range({domain1});

    EXPECT_EQ(range0, range);
    EXPECT_EQ(range1, range);

    Box domain = domain0 + domain1;

    auto X0 = forall_p<double, DIM>(f_iotaCorner, domain, dX);
    auto X_0 = forall<double, DIM>(f_shear, domain0, X0, stretch, angle);
    auto X_1 = forall<double, DIM>(f_shear, domain1, X0, stretch, angle);

    auto Adj0 = adj0(X_0);
    auto Adj1 = adj1(X_1);
    EXPECT_EQ(Adj0.box(), range);
    EXPECT_EQ(Adj1.box(), range);
    double Adj00 = Adj0.sum(0)/range.size();
    double Adj10 = Adj0.sum(1)/range.size();
    double Adj01 = Adj1.sum(0)/range.size();
    double Adj11 = Adj1.sum(1)/range.size();

    EXPECT_NEAR(Adj00, 1.0, 1e-12);
    EXPECT_NEAR(Adj10, -b, 1e-12);
    EXPECT_NEAR(Adj01, 0, 1e-12);
    EXPECT_NEAR(Adj11, a, 1e-12);

    EXPECT_NEAR(sqrt(Adj00*Adj00 + Adj10*Adj10), da[0], 1e-12);
    EXPECT_NEAR(sqrt(Adj01*Adj01 + Adj11*Adj11), da[1], 1e-12);

    auto Adj0_ = Operator::cofactor(X_0, 0);
    auto Adj1_ = Operator::cofactor(X_1, 1);

    Adj0_ -= Adj0;
    Adj1_ -= Adj1;

    EXPECT_LT(Adj0_.absMax(), 1e-12);
    EXPECT_LT(Adj1_.absMax(), 1e-12);

    #if PR_VERBOSE > 0
    HDF5Handler h5;
    h5.writePatch(X0, "TEST_ADJUGATE_X0");
    h5.writePatch(X_0, "TEST_ADJUGATE_X_0");
    h5.writePatch(X_1, "TEST_ADJUGATE_X_1");
    h5.writePatch(Adj0, "TEST_ADJUGATE_ADJ0");
    h5.writePatch(Adj1, "TEST_ADJUGATE_ADJ1");
    h5.writePatch(Adj0_, "TEST_ADJUGATE_ADJ0_");
    h5.writePatch(Adj1_, "TEST_ADJUGATE_ADJ1_");
    #endif
}
#endif

int main(int argc, char *argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
#ifdef PR_MPI
    MPI_Init(&argc, &argv);
#endif
    std::cout << "DIM = " << DIM << std::endl;
    int result = RUN_ALL_TESTS();
#ifdef PR_MPI
    MPI_Finalize();
#endif
    return result;
}
