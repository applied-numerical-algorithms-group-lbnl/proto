#include <gtest/gtest.h>
#include "Proto.H"
#include "Lambdas.H"

TEST(InterpStencil, Constant) {
    
    constexpr unsigned int C = 3;
    constexpr unsigned char D = 1;

    int domainSize = 8;
    Point refRatio(2,4,2,4,2,4);

    Box srcBox = Box::Cube(domainSize);
    Box dstBox = srcBox.refine(refRatio);

    BoxData<double, C, HOST, D> hostSrcData(srcBox);
    BoxData<double, C, HOST, D> hostDstData(dstBox);
    BoxData<double, C, HOST, D> hostSlnData(dstBox);
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

    I.apply(hostDstData, hostSrcData);
    EXPECT_EQ(I.span(), Box(Point::Ones()));
    EXPECT_EQ(I.ghost(), Point::Zeros());
    hostSlnData -= hostDstData;
    EXPECT_LT(hostSlnData.absMax(), 1e-12);
}

TEST(InterpStencil, Linear) {
    HDF5Handler h5;

    constexpr unsigned int C = 3;
    constexpr unsigned char D = 1;

    int domainSize = 128;
    Point refRatio(2,4,2,4,2,4);
    Array<double, DIM> k{1,2,3,4,5,6};
    Array<double, DIM> offset{1,2,3,4,5,6};

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
        I.apply(hostDstData, hostSrcData);
        
        hostSlnData -= hostDstData;
        err[nn] = hostSlnData.absMax();
        domainSize *= 2;
    }
    
    for (int ii = 1; ii < N; ii++)
    {
        double rate = log(err[ii-1]/err[ii]) / log(2.0);
        EXPECT_LT(std::abs(1-rate), 0.01);
    }
}

TEST(InterpStencil, Quadratic) {
    HDF5Handler h5;

    constexpr unsigned int C = 3;
    constexpr unsigned char D = 1;

    int domainSize = 128;
    Point refRatio(2,4,2,4,2,4);
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

        I.apply(hostDstData, hostSrcData);
        
        hostSlnData -= hostDstData;
        err[nn] = hostSlnData.absMax();
        PR_DEBUG_MSG(1, "error: %3.2e", err[nn]);
        domainSize *= 2;
    }
    
    for (int ii = 1; ii < N; ii++)
    {
        double rate = log(err[ii-1]/err[ii]) / log(2.0);
        PR_DEBUG_MSG(1,"rate: %3.2f", rate);
        EXPECT_GT(rate, 2 - 0.01);
    }
}

TEST(InterpStencil, FiniteVolume) {
    HDF5Handler h5;

    constexpr unsigned int C = 3;
    constexpr unsigned char D = 1;

    int domainSize = 32;
    Point refRatio(2,4,2,4,2,4);
    Array<double, DIM> k{1,2,3,4,5,6};
    Array<double, DIM> offset{1,2,3,4,5,6};
    offset *= (0.1);

    int N = 4;
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
        BoxData<double, C, HOST, D> hostDstData4(dstBox);
        BoxData<double, C, HOST, D> hostDstData5(dstBox);
        BoxData<double, C, HOST, D> hostSlnData(dstBox);
        BoxData<double, C, HOST, D> hostErrData4(dstBox);
        BoxData<double, C, HOST, D> hostErrData5(dstBox);
        forallInPlace_p(f_phi_avg, hostSrcData, cdx, k, offset);
        forallInPlace_p(f_phi_avg, hostSlnData, fdx, k, offset);

        auto I4 = InterpStencil<double>::FiniteVolume(refRatio, 4);
        auto I5 = InterpStencil<double>::FiniteVolume(refRatio, 5);
        EXPECT_EQ(I4.span(), Box::Kernel(2));
        EXPECT_EQ(I5.span(), Box::Kernel(2));
        EXPECT_EQ(I4.ghost(), Point::Ones(2));
        EXPECT_EQ(I5.ghost(), Point::Ones(2));

        I4.apply(hostDstData4, hostSrcData);
        I5.apply(hostDstData5, hostSrcData);
        
        hostDstData4.copyTo(hostErrData4);
        hostDstData5.copyTo(hostErrData5);
        hostErrData4 -= hostSlnData;
        hostErrData5 -= hostSlnData;
        err4[nn] = hostErrData4.absMax();
        err5[nn] = hostErrData4.absMax();
        PR_DEBUG_MSG(1, "4th Order Error: %3.2e", err4[nn]);
        PR_DEBUG_MSG(1, "5th Order Error: %3.2e", err5[nn]);
        domainSize *= 2;
    }
    
    for (int ii = 1; ii < N; ii++)
    {
        double rate4 = log(err4[ii-1]/err4[ii]) / log(2.0);
        double rate5 = log(err5[ii-1]/err5[ii]) / log(2.0);
        PR_DEBUG_MSG(1,"4th Order Convergence Rate: %3.2f", rate4);
        PR_DEBUG_MSG(1,"5th Order Convergence Rate: %3.2f", rate5);
        EXPECT_GT(rate4, 4 - 0.01);
        // Currently not working
        //EXPECT_GT(rate5, 5 - 0.01);
    }
}
#if 0
PROTO_KERNEL_START
void srcLambdaF(Point p, Var<double>& v, const double dx) {
    //v(0) = p[0]*dx;
    v(0) = sin(dx*p[0]);
    #if DIM > 1
    v(0) *= cos(dx*p[1]);
    #endif
}
PROTO_KERNEL_END(srcLambdaF, srcLambda)

PROTO_KERNEL_START
void solnLambdaF(Point p, Var<double>& v, const double dx) {
    //v(0) = p[0]*dx/3.0;
    v(0) = sin(dx/3.0*p[0]);
    #if DIM > 1
    v(0) *= cos(dx/3.0*p[1]);
    #endif
}
PROTO_KERNEL_END(solnLambdaF, solnLambda)

TEST(Interp, PiecewiseLinear) {
    Point r = Point::Ones(3);
    auto PWC = InterpStencil<double>::PiecewiseConstant(r);
    auto PWL = InterpStencil<double>::PiecewiseLinear(r);

    int domainSize = 16, numIter = 4;
    double error_C[numIter];
    double error_L[numIter];
    for (int ii = 0; ii < numIter; ii++) {
        Box B0 = Box::Cube(domainSize);
        Box B1 = Box(B0.low()*r, B0.high()*r);
        Box B2 = B0.refine(r);
        BoxData<double> Src(B0);
        BoxData<double> DC0(B2,1337.);
        BoxData<double> DL0(B1,1337.);
        BoxData<double> DC1(B2,17.);
        BoxData<double> DL1(B1,17.);
        BoxData<double> Soln(B2);

        //double dx = (M_PI/4.0)/domainSize;
        double dx = 1.0/domainSize;


        forallInPlace_p(srcLambda,Src,dx);
        forallInPlace_p(solnLambda,Soln,dx);

        DC0 |= PWC(Src);
        DL0 |= PWL(Src);
        DC1 += PWC(Src);
        DL1 += PWL(Src);
        BoxData<double> DC2 = PWC(Src);
        BoxData<double> DL2 = PWL(Src);

        EXPECT_EQ(DC2.box(),B2);
        EXPECT_EQ(DL2.box(),B1);
        DC0 -= DC2;
        DL0 -= DL2;
        DC1 -= 17;
        DL1 -= 17;
        DC1 -= DC2;
        DL1 -= DL2;

        EXPECT_DOUBLE_EQ(DC0.absMax(),0.);
        EXPECT_DOUBLE_EQ(DL0.absMax(),0.);
        EXPECT_LT(DC1.absMax(),1e-6);
        EXPECT_LT(DL1.absMax(),1e-6);
        DC2 -= Soln;
        DL2 -= Soln;
        error_C[ii] = DC2.absMax();
        error_L[ii] = DL2.absMax();
        domainSize *= 2;
    }

    for (int ii = 1; ii < numIter; ii++) {
        double rate = log2(error_C[ii-1]/error_C[ii]);
        EXPECT_LT(abs(rate - 1.),0.1);
        rate = log2(error_L[ii-1]/error_L[ii]);
        EXPECT_LT(abs(rate - 2.),0.1);
    }
} 

TEST(InterpStencil, BoxInference) {
    Box B0 = Box::Cube(4).grow(1);
    Box B1 = Box::Cube(8).grow(1);
    Box K = Box(Point::Ones(-2),Point::Ones(8));
    auto Src = forall_p<double>(pointSum, B0);
    auto Soln = forall_p<double>(halfPointSum, K);

    BoxData<double> Dest0(B1,1337.);
    BoxData<double> Dest1(B1,1337.);
    BoxData<double> Dest2(B1,1337.);
    auto I = InterpStencil<double>::PiecewiseLinear(Point::Ones(2));

    Dest0 |= I(Src, B0.grow(-1));
    Dest1 |= I(Src, B0.grow(1));
    Dest2 |= I(Src);
    BoxData<double> Dest3 = I(Src);

    EXPECT_EQ(Dest3.box(),K);

    BoxData<double,1,HOST> D0h(B1), D1h(B1), D2h(B1), D3h(K), Solnh(K);
    Dest0.copyTo(D0h); Dest1.copyTo(D1h); Dest2.copyTo(D2h); Dest3.copyTo(D3h);
    Soln.copyTo(Solnh);

    for (auto iter : K) {
        EXPECT_EQ(D3h(iter),Solnh(iter));
        if (B1.grow(-1).contains(iter))
            EXPECT_EQ(D0h(iter),Solnh(iter));
        if (B1.contains(iter)) {
            EXPECT_EQ(D1h(iter),Solnh(iter));
            EXPECT_EQ(D2h(iter),Solnh(iter));
        }
    }
}

TEST(InterpStencil, CosApply) {
    int domainSize = 8, numIter = 5;
    for (int ii = 0; ii < numIter; ii++) {
        Box B0 = Box::Cube(domainSize).grow(1);
        Box B1 = Box::Cube(2*domainSize).grow(1);
        double dx = M_PI/domainSize;
        BoxData<double> Src  = forall_p<double>(cosxCosyFunc,     B0, dx);
        BoxData<double> Soln = forall_p<double>(cosxCosyPCosFunc, B1, dx);
        BoxData<double> Dest = forall_p<double>(cosxFunc,         B1, dx);

        auto interp = InterpStencil<double>::PiecewiseLinear(Point::Ones(2));
        Dest += interp(Src);
        BoxData<double,1,HOST> dest_host(B1), src_host(B0);
        Dest.copyTo(dest_host); Src.copyTo(src_host); 
        for (auto iter : B1) {
            double x =iter[0]*dx/2.0;
            double y =iter[1]*dx/2.0;
            Point p = iter % Point::Ones(2);
            Point q = iter / Point::Ones(2);
            double value = cos(x);
            if ( p == Point::Zeros()) {
                value += src_host(q);
                EXPECT_FLOAT_EQ(dest_host(iter),value);
            } else if ((p == Point::Basis(0)) ||
                       (p == Point::Basis(1)) ||
                       (p == Point::Basis(0,-1)) ||
                       (p == Point::Basis(1,-1))) {
                value += (src_host(q) + src_host(q+p))/2.0;
                EXPECT_FLOAT_EQ(dest_host(iter),value);
            } else if (p ==  Point::Ones()) {
                value += (src_host(q) + src_host(q+p) + src_host(q + Point::Basis(0)) + src_host(q + Point::Basis(1)))/4.0;
                EXPECT_FLOAT_EQ(dest_host(iter),value);
            }
        }
        Dest -= Soln;
        domainSize *= 2;
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
