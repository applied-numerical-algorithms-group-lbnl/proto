#include <gtest/gtest.h>
#include "Proto.H"
#include "TestFunctions.H"
#include "BoxOp_TestLaplace.H"

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

TEST(LevelOp, LaplaceApply) {
    HDF5Handler h5;

    int domainSize = 64;
    double a0 = 0.125;
    Array<double, DIM> k{1,1,1,1,1,1};
    Array<double, DIM> a{a0, a0, a0, a0, a0, a0};
    typedef BoxOp_TestLaplace<double> OP;
    
    int N = 3;
    double err[N];
    for (int nn = 0; nn < N; nn++)
    {
        Point boxSize = Point::Ones(domainSize / 4);
        Array<double, DIM> dx = Point::Ones();
        dx /= domainSize;
        
        auto layout = testLayout(domainSize, boxSize);

        LevelOp<BoxOp_TestLaplace, double> op(layout, dx);
        
        LevelBoxData<double, 1> srcData(layout, OP::ghost());
        LevelBoxData<double, 1> dstData(layout, Point::Zeros());
        LevelBoxData<double, 1> slnData(layout, Point::Zeros());
        LevelBoxData<double, 1> errData(layout, Point::Zeros());
        srcData.initialize(f_phi,  dx, k, a);
        slnData.initialize(f_Lphi, dx, k, a);
        dstData.setVal(0);
        errData.setVal(0);

        op(dstData, srcData);

        dstData.copyTo(errData);
        errData.increment(slnData, -1.0);
        err[nn] = errData.absMax(); 
#if PR_VERBOSE > 0
        std::cout << "Error (Max Norm): " << err[nn] << std::endl;
        h5.writeLevel({"phi"},  dx, srcData, "TEST_LEVELOP_SRC_N%i", nn);
        h5.writeLevel({"Lphi"}, dx, dstData, "TEST_LEVELOP_DST_N%i", nn);
        h5.writeLevel({"Lphi"}, dx, slnData, "TEST_LEVELOP_SLN_N%i", nn);
        h5.writeLevel({"error"}, dx, errData, "TEST_LEVELOP_ERR_N%i", nn);
#endif
        domainSize *= 2;
    }

    for (int ii = 1; ii < N; ii++)
    {
        double rate = log(err[ii-1]/err[ii])/log(2.0);
        double rateErr = std::abs(rate - 4.0);
#if PR_VERBOSE > 0
        std::cout << "Convergence rate: " << rate << std::endl;
#endif
        EXPECT_LT(rateErr, 0.3);
    }
}

TEST(LevelOp, ApplyBCPeriodic) {
    HDF5Handler h5;

    int domainSize = 64;
    double a0 = 0.125;
    Array<double, DIM> k{1,1,1,1,1,1};
    Array<double, DIM> a{a0, a0, a0, a0, a0, a0};
    typedef BoxOp_TestLaplace<double> OP;
   
    Box domainBox = Box::Cube(domainSize); 
    std::array<bool, DIM> periodicity;
    periodicity.fill(true);
    Point boxSize = Point::Ones(domainSize / 4);
    ProblemDomain domain(domainBox, periodicity);
    DisjointBoxLayout layout(domain, boxSize);

    Array<double, DIM> dx = Point::Ones();
    dx /= domainSize;
    LevelOp<BoxOp_TestLaplace, double> op(layout, dx);
    
    LevelBoxData<double, 1> srcData0(layout, Point::Zeros());
    LevelBoxData<double, 1> srcData1(layout, OP::ghost());
    LevelBoxData<double, 1> srcData2(layout, OP::ghost());
    LevelBoxData<double, 1> errData(layout, OP::ghost());
    
    srcData0.initialize(f_phi, dx, k, a);
    srcData2.initialize(f_phi, dx, k, a);

    srcData1.setVal(0.0);
    for (auto iter : layout)
    {
        srcData0[iter].copyTo(srcData1[iter]);
    }
#if PR_VERBOSE > 0
    h5.writeLevel({"phi"},  dx, srcData0, "TEST_LEVELOP_BC_PHI_0");
    h5.writeLevel({"phi"},  dx, srcData1, "TEST_LEVELOP_BC_PHI_1");
#endif
    op.applyBC(srcData1);
    errData.setVal(0);

    for (auto iter : layout)
    {
        srcData1[iter].copyTo(errData[iter]);
        errData[iter] -= srcData2[iter];
    }
    double error = errData.absMax();
#if PR_VERBOSE > 0
    std::cout << "Error (Max Norm): " << error << std::endl;
    h5.writeLevel({"phi"},  dx, srcData1, "TEST_LEVELOP_BC_PHI_2");
    h5.writeLevel({"err"},  dx, errData,  "TEST_LEVELOP_BC_ERR");
#endif
    EXPECT_LT(error, 1e-12);
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
