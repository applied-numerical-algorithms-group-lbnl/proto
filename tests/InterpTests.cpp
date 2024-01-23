#include <gtest/gtest.h>
#include "Proto.H"
#include "Lambdas.H"

TEST(InterpStencil, BoxInference)
{
    int domainSize = 8;
    Point refRatio(2,4,2,4,2,4);
    auto I0 = InterpStencil<double>::Constant(refRatio);
    auto I1 = InterpStencil<double>::Linear(refRatio);
    auto I2 = InterpStencil<double>::Quadratic(refRatio);

    Box B = Box::Cube(domainSize);
    Box B0 = B.refine(refRatio);
    Box B1 = B.extrude(Point::Ones(), -1).refine(refRatio);
    Box B2 = B.grow(-1).refine(refRatio);
    
    Box b0 = B0.grow(-1);
    Box b1 = B1.grow(-1);
    Box b2 = B2.grow(-1);

    BoxData<double> hostSrcData(B);
    
    BoxData<double> hostDstData_0 = I0(hostSrcData);
    BoxData<double> hostDstData_1 = I1(hostSrcData);
    BoxData<double> hostDstData_2 = I2(hostSrcData);
    
    BoxData<double> hostDstData_3 = I0(hostSrcData, b0);
    BoxData<double> hostDstData_4 = I1(hostSrcData, b1);
    BoxData<double> hostDstData_5 = I2(hostSrcData, b2);

    EXPECT_EQ(hostDstData_0.box(), B0);
    EXPECT_EQ(hostDstData_1.box(), B1);
    EXPECT_EQ(hostDstData_2.box(), B2);

    EXPECT_EQ(hostDstData_3.box(), b0);
    EXPECT_EQ(hostDstData_4.box(), b1);
    EXPECT_EQ(hostDstData_5.box(), b2);

}
TEST(InterpStencil, Constant) {
    
    constexpr unsigned int C = 3;
    constexpr unsigned int D = 1;

    int domainSize = 2;
    Point refRatio(2,4,2,4,2,4);

    Box srcBox = Box::Cube(domainSize);
    Box dstBox = srcBox.refine(refRatio);

    BoxData<double, C, HOST, D> hostSrcData(srcBox);
    BoxData<double, C, HOST, D> hostDstData(dstBox);
    BoxData<double, C, HOST, D> hostSlnData(dstBox);
    BoxData<double, C, HOST, D> hostErrData(dstBox);
    
    forallInPlace_p(f_pointID, hostSrcData);
    
    for (auto p : srcBox)
    {
        Box b = Box(p,p).refine(refRatio);
        for (auto q : b)
        {
            for (int cc = 0; cc < C; cc++)
            for (int dd = 0; dd < D; dd++)
            {
                hostSlnData(q,cc,dd) = hostSrcData(p,cc,dd);
            }
        }
    }
    
    auto I = InterpStencil<double>::Constant(refRatio);
    EXPECT_EQ(I.span(), Box(Point::Ones()));
    EXPECT_EQ(I.ghost(), Point::Zeros());
    
    hostDstData |= I(hostSrcData);
    hostDstData.copyTo(hostErrData);
    hostErrData -= hostSlnData;
    EXPECT_LT(hostErrData.absMax(), 1e-12);
    
    BoxData<double, C, HOST, D> hostOutData = I(hostSrcData);
    EXPECT_EQ(hostOutData.box(), dstBox);
    hostOutData.copyTo(hostErrData);
    hostErrData -= hostSlnData;
    EXPECT_LT(hostErrData.absMax(), 1e-12);
    
    Box limitBox = dstBox.grow(-Point::Basis(0));
    hostDstData.setVal(7);
    hostDstData.printData();
    hostDstData |= I(hostSrcData, limitBox);
    for (auto pi : hostDstData.box())
    {
        for (int cc = 0; cc < C; cc++)
        for (int dd = 0; dd < D; dd++)
        {
            if (limitBox.contains(pi))
            {
                EXPECT_EQ(hostDstData(pi,cc,dd), hostSlnData(pi, cc, dd));
            } else {
                EXPECT_EQ(hostDstData(pi,cc,dd), 7);
            }
        }
    }
    
    hostDstData.setVal(7);
    hostDstData += I(hostSrcData);
    hostSlnData += 7;
    hostDstData.copyTo(hostErrData);
    hostErrData -= hostSlnData;
    EXPECT_LT(hostErrData.absMax(), 1e-12);
    
}
#if 0
TEST(InterpStencil, Linear) {
#ifdef PR_HDF5
	HDF5Handler h5;
#endif

    constexpr unsigned int C = 3;
    constexpr unsigned int D = 1;

    int domainSize = 128;
    Point refRatio(2,4,2,4,2,4);
    Array<double, DIM> k{1,2,3,4,5,6};
    Array<double, DIM> offset{1,2,3,4,5,6};

    auto I0 = InterpStencil<double>::Linear(refRatio);
    Box b0 = Box::Cube(8);
    Box domainSln = b0.grow(PR_NODE);
    Box range = I0.range(domainSln);
    Box rangeSln = b0.refine(refRatio);
    Box domain = I0.domain(rangeSln); 
    EXPECT_EQ(range,rangeSln);
    EXPECT_EQ(domain,domainSln);

    int N = 2;
    double err[N];
    for (int nn = 0; nn < N; nn++)
    {
        Array<double, DIM> cdx(1.0/domainSize);
        Array<double, DIM> fdx = cdx;
        fdx /= (Array<double, DIM>)refRatio;
        Box srcBox = Box::Cube(domainSize);
        Box dstBox = srcBox.refine(refRatio);

        BoxData<double, C, HOST, D> hostSrcData(srcBox.grow(PR_NODE));
        BoxData<double, C, HOST, D> hostDstData(dstBox);
        BoxData<double, C, HOST, D> hostSlnData(dstBox);
        forallInPlace_p(f_phi, hostSrcData, cdx, k, offset);
        forallInPlace_p(f_phi, hostSlnData, fdx, k, offset);

        auto I = InterpStencil<double>::Linear(refRatio);
        EXPECT_EQ(I.span(), Box(Point::Ones(2)));
        EXPECT_EQ(I.ghost(), Point::Ones());
        hostDstData |= I(hostSrcData);
        
        hostSlnData -= hostDstData;
        err[nn] = hostSlnData.absMax();
        domainSize *= 2;
    }
    
    for (int ii = 1; ii < N; ii++)
    {
        double rate = log(err[ii-1]/err[ii]) / log(2.0);
        EXPECT_GT(rate, 2.0 - 0.01);
    }
}
#endif
#if 0
TEST(InterpStencil, Quadratic) {
#ifdef PR_HDF5
    HDF5Handler h5;
#endif

    constexpr unsigned int C = 3;
    constexpr unsigned int D = 1;

    int domainSize = 128;
    Point refRatio(2,4,2,4,2,4);
    
    auto I0 = InterpStencil<double>::Quadratic(refRatio);

    Box b0 = Box::Cube(8);
    Box domainSln = b0.grow(1);
    Box range = I0.range(domainSln);
    Box rangeSln = b0.refine(refRatio);
    Box domain = I0.domain(rangeSln); 
    EXPECT_EQ(range,rangeSln);
    EXPECT_EQ(domain,domainSln);
    
    Array<double, DIM> k{1,2,3,4,5,6};
    Array<double, DIM> offset{1,2,3,4,5,6};
    offset *= (0.1);
    int N = 2;
    double err[N];
    for (int nn = 0; nn < N; nn++)
    {
        Array<double, DIM> cdx(1.0/domainSize);
        Array<double, DIM> fdx = cdx;
        fdx /= (Array<double, DIM>)refRatio;
        Box srcBox = Box::Cube(domainSize);
        Box dstBox = srcBox.refine(refRatio);

        BoxData<double, C, HOST, D> hostSrcData(srcBox.grow(1));
        BoxData<double, C, HOST, D> hostDstData(dstBox);
        BoxData<double, C, HOST, D> hostSlnData(dstBox);
        forallInPlace_p(f_phi, hostSrcData, cdx, k, offset);
        forallInPlace_p(f_phi, hostSlnData, fdx, k, offset);

        auto I = InterpStencil<double>::Quadratic(refRatio);
        EXPECT_EQ(I.span(), Box::Kernel(1));
        EXPECT_EQ(I.ghost(), Point::Ones());

        hostDstData |= I(hostSrcData);
        
        hostSlnData -= hostDstData;
        err[nn] = hostSlnData.absMax();
        PR_DEBUG_MSG(1, "Error: %3.2e", err[nn]);
        domainSize *= 2;
    }
    
    for (int ii = 1; ii < N; ii++)
    {
        double rate = log(err[ii-1]/err[ii]) / log(2.0);
        PR_DEBUG_MSG(1,"Rate: %3.2f", rate);
        EXPECT_GT(rate, 3 - 0.01);
    }
}
#endif
TEST(InterpStencil, FiniteVolume) {
#ifdef PR_HDF5
    HDF5Handler h5;
#endif

    constexpr unsigned int C = 3;
    constexpr unsigned int D = 1;

    int domainSize = 128;
    Point refRatio(2,4,2,4,2,4);
    Array<double, DIM> k{1,2,3,4,5,6};
    Array<double, DIM> offset{1,2,3,4,5,6};
    offset *= (0.1);

    int N = 2;
    double err2[N];
    double err3[N];
    double err4[N];
    double err5[N];
    for (int nn = 0; nn < N; nn++)
    {
        Array<double, DIM> cdx(1.0/domainSize);
        Array<double, DIM> fdx = cdx;
        fdx /= (Array<double, DIM>)refRatio;
        Box srcBox = Box::Cube(domainSize);
        Box dstBox = srcBox.refine(refRatio);

        BoxData<double, C, HOST, D> hostSrcData(srcBox.grow(2));
        BoxData<double, C, HOST, D> hostDstData2(dstBox);
        BoxData<double, C, HOST, D> hostDstData3(dstBox);
        BoxData<double, C, HOST, D> hostDstData4(dstBox);
        BoxData<double, C, HOST, D> hostDstData5(dstBox);
        BoxData<double, C, HOST, D> hostSlnData(dstBox);
        BoxData<double, C, HOST, D> hostErrData2(dstBox);
        BoxData<double, C, HOST, D> hostErrData3(dstBox);
        BoxData<double, C, HOST, D> hostErrData4(dstBox);
        BoxData<double, C, HOST, D> hostErrData5(dstBox);
        forallInPlace_p(f_phi_avg, hostSrcData, cdx, k, offset);
        forallInPlace_p(f_phi_avg, hostSlnData, fdx, k, offset);
        auto I2 = InterpStencil<double>::FiniteVolume(refRatio, 2);
        auto I3 = InterpStencil<double>::FiniteVolume(refRatio, 3);
        auto I4 = InterpStencil<double>::FiniteVolume(refRatio, 4);
        auto I5 = InterpStencil<double>::FiniteVolume(refRatio, 5);
        EXPECT_EQ(I2.span(), Box::Kernel(2));
        EXPECT_EQ(I3.span(), Box::Kernel(2));
        EXPECT_EQ(I2.ghost(), Point::Ones(2));
        EXPECT_EQ(I3.ghost(), Point::Ones(2));
        EXPECT_EQ(I4.span(), Box::Kernel(2));
        EXPECT_EQ(I5.span(), Box::Kernel(2));
        EXPECT_EQ(I4.ghost(), Point::Ones(2));
        EXPECT_EQ(I5.ghost(), Point::Ones(2));

        hostDstData2 |= I2(hostSrcData);
        hostDstData3 |= I3(hostSrcData);
        hostDstData4 |= I4(hostSrcData);
        hostDstData5 |= I5(hostSrcData);

        hostDstData2.copyTo(hostErrData2);
        hostDstData3.copyTo(hostErrData3);
        hostDstData4.copyTo(hostErrData4);
        hostDstData5.copyTo(hostErrData5);
        hostErrData2 -= hostSlnData;
        hostErrData3 -= hostSlnData;
        err2[nn] = hostErrData2.absMax();
        err3[nn] = hostErrData3.absMax();
        PR_DEBUG_MSG(1, "2nd Order Error: %3.2e", err2[nn]);
        PR_DEBUG_MSG(1, "3rd Order Error: %3.2e", err3[nn]);
        hostErrData4 -= hostSlnData;
        hostErrData5 -= hostSlnData;
        err4[nn] = hostErrData4.absMax();
        err5[nn] = hostErrData5.absMax();
        PR_DEBUG_MSG(1, "4th Order Error: %3.2e", err4[nn]);
        PR_DEBUG_MSG(1, "5th Order Error: %3.2e", err5[nn]);
        domainSize *= 2;
    }
    
    for (int ii = 1; ii < N; ii++)
    {
        double rate2 = log(err2[ii-1]/err2[ii]) / log(2.0);
        double rate3 = log(err3[ii-1]/err3[ii]) / log(2.0);
        PR_DEBUG_MSG(1,"2nd Order Convergence Rate: %3.2f", rate2);
        PR_DEBUG_MSG(1,"3rd Order Convergence Rate: %3.2f", rate3);
        EXPECT_GT(rate2, 2 - 0.01);
        EXPECT_GT(rate3, 3 - 0.01);
        double rate4 = log(err4[ii-1]/err4[ii]) / log(2.0);
        double rate5 = log(err5[ii-1]/err5[ii]) / log(2.0);
        PR_DEBUG_MSG(1,"4th Order Convergence Rate: %3.2f", rate4);
        PR_DEBUG_MSG(1,"5th Order Convergence Rate: %3.2f", rate5);
        EXPECT_GT(rate4, 4 - 0.01);
        EXPECT_GT(rate5, 5 - 0.01);
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
