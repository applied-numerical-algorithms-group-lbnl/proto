#include <gtest/gtest.h>
#include "Proto.H"
#include "Lambdas.H"
#include <iomanip>

#define NCOMP 1

using namespace Proto;
template<typename T, MemType MEM>
PROTO_KERNEL_START
void f_computeExp_tmp(Var<T, 1, MEM>& a_xp, const Var<T, DIM, MEM>& a_x, Point a_p)
{
    a_xp(0) = 1.0;
    for (int dir = 0; dir < DIM; dir++)
    {
        a_xp(0) *= pow(a_x(dir), a_p[dir]);
    }
}
PROTO_KERNEL_END(f_computeExp_tmp, f_computeExp)

template<typename MAP>
std::vector<Matrix<double>> computeM(MBMap<MAP>& a_map, int a_boxSize, Point a_dst)
{
    double h = 1.0/a_boxSize;
    auto C2C = Stencil<double>::CornersToCells(4);
    
    Box Bx = Box::Cube(a_boxSize*2);
    Box Bx0 = C2C.domain(Bx).extrude(Point::Ones(),-1);
   
    std::vector<BoxData<double, DIM>> X;
    std::vector<Box> B;
    auto Xg = a_map(a_dst, Point::Zeros(), 0);
    for (int ii = 0; ii < 4; ii++)
    {
        X.push_back(a_map(Bx0,ii,0));
        X[ii] -= Xg;
    }
    B.push_back(Box::Cube(a_boxSize));
    B.push_back(B[0].shift(Point::X()*a_boxSize));
    B.push_back(B[1].shift(Point::Y()*a_boxSize));
    B.push_back(B[0].shift(Point::Y()*a_boxSize));

    std::vector<Point> exponents;
    for (auto pi : Box::Cube(3))
    {
        if (pi.sum() <= 2)
        {
            exponents.push_back(pi);
        }
    }
   
    std::vector<std::vector<BoxData<double>>> Xp;
    Xp.resize(4);
    std::vector<std::vector<BoxData<double>>> Xp_avg;
    Xp_avg.resize(4);
    for (int bi = 0; bi < 4; bi++)
    {
        Xp[bi].resize(exponents.size());
        for (int ei = 0; ei < exponents.size(); ei++)
        {
            Xp[bi][ei].define(X[bi].box());
            forallInPlace(f_computeExp, Xp[bi][ei], X[bi], exponents[ei]);
            Xp_avg[bi].push_back(C2C(Xp[bi][ei]));
        }
    }

    std::vector<Point> footprint;
    for (auto pi : Box::Kernel(2))
    {
        if (pi.abs().sum() <= 2)
        {
            footprint.push_back(pi + a_dst);
        }
    }
    std::vector<std::tuple<Point, int>> srcs;
    for (int si = 0; si < footprint.size(); si++)
    {
        Point s = footprint[si];
        int block = -1;
        for (int bi = 0; bi < 4; bi++)
        {
            if (B[bi].contains(s)){ block = bi; break; }
        }
        if (block >= 0)
        {
            srcs.push_back(std::make_tuple(s, block));
        }
        
    }

    std::vector<Matrix<double>> M(2);
    M[0].define(srcs.size(), exponents.size());
    for (int si = 0; si < srcs.size(); si++)
    {
        Point s = std::get<0>(srcs[si]);
        int block = std::get<1>(srcs[si]);
        auto& src = srcs[si];
        for (auto ei = 0; ei < exponents.size(); ei++)
        {
            M[0](si, ei) = Xp_avg[block][ei](s);
            //std::cout << std::setw(10) << std::setprecision(5) << Xp_avg[block][ei](s) << ", ";
        }
        //std::cout << std::endl;
    }
 
    M[1].define(1, exponents.size());
    //std::cout << "\n" << std::endl;
    for (auto ei = 0; ei < exponents.size(); ei++)
    {
        M[1](0, ei) = Xp_avg[0][ei](a_dst);
        //std::cout << std::setw(10) << std::setprecision(5) << Xp_avg[0][ei](a_dst) << ", ";
    }
    //std::cout << std::endl;

    return M;
}

#if 0
TEST(MBPointInterpOp, CheckMatrix)
{

    auto C2C = Stencil<double>::CornersToCells(4);

    int domainSize = 8;
    int boxSize = 8;
    
    auto domain = buildShear(domainSize);
    Point boxSizeVect = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);

    Array<Point, DIM+1> ghost;
    ghost.fill(Point::Ones(4));
    ghost[0] = Point::Ones(1);
    Point boundGhost = Point::Ones();

    // initialize map
    MBMap<ShearMap_t> map(ShearMap, layout, ghost, boundGhost);

    // input footprint
    std::vector<Point> footprint;
    for (auto pi : Box::Kernel(1))
    {
        footprint.push_back(pi);
        if (pi.codim() == 1)
        {
            footprint.push_back(pi*2);
        }
    }
    
    MBInterpOp interp(map, footprint, ghost[0], 4);
    
    auto M = computeM(map, boxSize, Point::X()*boxSize);
    Box B0 = Box::Cube(boxSize);
    double CErr = 0;
    double DErr = 0;
    double SErr = 0;
    for (auto pi : B0.extrude(Point::Ones(), 1))
    {
        if (B0.contains(pi)) { continue; }
        std::cout << "Checking point " << pi << std::endl;
        MBDataPoint t(*layout.begin(), pi, layout);
        auto& op = interp(t);
        auto M = computeM(map, boxSize, pi);
        auto Cinv = M[0].inverse();
        auto S = M[1]*Cinv;
        auto EC = M[0] - op.MC();
        auto ED = M[1] - op.MD();
        auto ES = S - op.MS();
        CErr = std::max(CErr, (EC.absMax()));
        DErr = std::max(DErr, (ED.absMax()));
        SErr = std::max(SErr, (ES.absMax()));
        M[0].print();
        M[1].print();
        S.print();
    }
    std::cout << "Error in C: " << CErr << std::endl;
    std::cout << "Error in D: " << DErr << std::endl;
    std::cout << "Error in S: " << SErr << std::endl;
}
#endif
#if 1
TEST(MBPointInterpOp, ShearApply) {
    int domainSize = 8;
    int boxSize = 8;
    HDF5Handler h5;

    constexpr int N = 3;
    double errInf[N];
    double errL1[N];
    for (int nn = 0; nn < N; nn++)
    {
        errInf[nn] = 0.0;
        errL1[nn] = 0.0;
        auto domain = buildShear(domainSize);
        Point boxSizeVect = Point::Ones(boxSize);
        MBDisjointBoxLayout layout(domain, boxSizeVect);

        Array<Point, DIM+1> ghost;
        ghost.fill(Point::Ones(4));
        ghost[0] = Point::Ones(1);
        Point boundGhost = Point::Ones();

        // initialize map
        MBMap<ShearMap_t> map(ShearMap, layout, ghost, boundGhost);

        // initialize data
        MBLevelBoxData<double, NCOMP, HOST> hostSrc(layout, ghost);
        MBLevelBoxData<double, NCOMP, HOST> hostDst(layout, ghost);
        MBLevelBoxData<double, NCOMP, HOST> hostErr(layout, ghost);
        MBLevelBoxData<double, 6, HOST> hostCoefs(layout, ghost);

        Array<double, DIM> exp{4,0,0,0,0,0};
        Array<double, DIM> offset{0,0,0,0,0,0};
        hostSrc.initConvolve(f_polyM, map, exp, offset);
        hostSrc.fillBoundaries();
        hostDst.setVal(0);
        hostErr.setVal(0);
        hostCoefs.setVal(0);
        
        // input footprint
        std::vector<Point> footprint;
        for (auto pi : Box::Kernel(1))
        {
            footprint.push_back(pi);
        }
        for (int dir = 0; dir < DIM; dir++)
        {
            footprint.push_back(Point::Basis(dir,2));
            footprint.push_back(Point::Basis(dir,-2));
        }

        // Create and Apply Operator
        MBInterpOp interp(map, footprint, ghost[0], 4);
        interp.apply(hostDst, hostSrc);
        interp.coefs(hostCoefs, hostSrc);
        Box domainBox = Box::Cube(domainSize);
        for (auto iter : layout)
        {
            auto block = layout.block(iter);
            Box B0 = layout[iter];
            auto& patch = hostDst[iter];
            auto& errPatch = hostErr[iter];
            for (auto dir : Box::Kernel(1))
            {
                if (!layout.domain().graph().isBlockBoundary(block, dir)) {continue; } 
                Box boundBox = B0.adjacent(dir*ghost[0]);
                if (domainBox.contains(boundBox)) {continue;}
                auto boundX = map(boundBox.grow(1), block, PR_CELL);
                BoxData<double, NCOMP, HOST> boundData0(boundBox.grow(1));
                BoxData<double, NCOMP, HOST> boundData(boundBox);
                forallInPlace_p(f_polyM, boundData0, block, boundX, exp, offset);
                Operator::convolve(boundData, boundData0);
                BoxData<double, NCOMP, HOST> boundErr(boundBox);
                patch.copyTo(boundErr);
                boundErr -= boundData;
                boundErr.copyTo(errPatch);
                double e = boundErr.absMax();
                errInf[nn] = e > errInf[nn] ? e : errInf[nn];
                errL1[nn] += boundErr.sum();
            }
        }
#if PR_VERBOSE > 0

        std::vector<std::string> coefNames;
        for (auto e : interp.exponents())
        {
            char var[100];
            sprintf(var, "%i_%i", e[0], e[1]);
            coefNames.push_back(var);
        }

        std::cout << "1/Jacobian: " << domainSize*domainSize << std::endl;
        h5.writeMBLevel(coefNames, map, hostCoefs, "MBInterpOpTests_Coefs_%i", nn);
        h5.writeMBLevel({"J"}, map, map.jacobian(), "MBInterpOpTests_J_%i", nn);
        h5.writeMBLevel({"phi"}, map, hostSrc, "MBInterpOpTests_Src_%i", nn);
        h5.writeMBLevel({"phi"}, map, hostDst, "MBInterpOpTests_Dst_%i", nn);
        h5.writeMBLevel({"err"}, map, hostErr, "MBInterpOpTests_Err_%i", nn);
        for (auto iter : layout)
        {
            pout() << "Block " << iter.global() << std::endl;
            hostSrc[iter].printData();
            hostDst[iter].printData(); 
        }
#endif
        PR_DEBUG_MSG(1, "Error (Max Norm): %3.2e", errInf[nn]);
        PR_DEBUG_MSG(1, "Error (L1 Norm): %3.2e", errL1[nn]);
        domainSize *= 2;
        boxSize *= 2;
    }

    double rateTol = 0.1;
    for (int ii = 1; ii < N; ii++)
    {
        double rateInf = log(errInf[ii-1]/errInf[ii])/log(2.0);
        double rateL1 = log(errL1[ii-1]/errL1[ii])/log(2.0);
        PR_DEBUG_MSG(1,"Rate (Max Norm): %3.2f", rateInf);
        PR_DEBUG_MSG(1,"Rate (L1 Norm): %3.2f", rateL1);
        double rateErr = abs(rateL1 - 4);
        EXPECT_LT(rateErr, rateTol);
    }
}
#endif
#if 0
TEST(MBPointInterpOp, XPointApply) {
    int domainSize = 16;
    int boxSize = 16;
    HDF5Handler h5;
    if (DIM > 2) 
    {
       // domainSize = 16;
       // boxSize = 8;
    }

    constexpr int N = 2;
    Array<double, N> err(0);
    for (int nn = 0; nn < N; nn++)
    {
        auto domain = buildXPoint(domainSize);
        Point boxSizeVect = Point::Ones(boxSize);
        MBDisjointBoxLayout layout(domain, boxSizeVect);

        Array<Point, DIM+1> ghost;
        ghost.fill(Point::Ones(3));
        ghost[0] = Point::Ones(1);
        Point boundGhost = Point::Ones();

        // initialize map
        MBMap<XPointMapRigid_t> map(XPointMapRigid, layout, ghost, boundGhost);

        // initialize data
        MBLevelBoxData<double, NCOMP, HOST> hostSrc(layout, ghost, boundGhost);
        MBLevelBoxData<double, NCOMP, HOST> hostDst(layout, ghost, boundGhost);
        MBLevelBoxData<double, NCOMP, HOST> hostErr(layout, ghost, boundGhost);

        Point k{1,1,1,1,1,1};
        Array<double, DIM> offset{1,1,1,1,1,1};
        hostSrc.initialize(f_phiM, map, k, offset);
        hostSrc.fillBoundaries();
        hostDst.setVal(0);
        hostErr.setVal(0);
#if PR_VERBOSE > 1
        h5.writeMBLevelBounds({"J"}, map.jacobian(), "MBInterpOpTests_J_N%i", nn);
#endif
        // input footprint
        std::vector<Point> footprint;
        for (auto pi : Box::Kernel(1))
        {
            footprint.push_back(pi);
        }
        for (int dir = 0; dir < DIM; dir++)
        {
            footprint.push_back(Point::Basis(dir,2));
            footprint.push_back(Point::Basis(dir,-2));
        }

        // Create and Apply Operator
        MBInterpOp interp(map, footprint, ghost[0], 4);
        interp.apply(hostDst, hostSrc);

        Box domainBox = Box::Cube(domainSize);
        for (auto iter : layout)
        {
            auto block = layout.block(iter);
            Box B0 = layout[iter];
            auto& patch = hostDst[iter];
            auto& errPatch = hostErr[iter];
            for (auto dir : Box::Kernel(1))
            {
                if (!layout.domain().graph().isBlockBoundary(block, dir)) {continue; } 
                //if (dir != Point::Basis(0) && dir != Point::Basis(1)
                //        && dir != (Point::Basis(1) + Point::Basis(0))){continue;}
                Box boundBox = B0.adjacent(dir*ghost[0]);
                if (domainBox.contains(boundBox)) {continue;}
                auto boundX = map(boundBox, block, PR_CELL);
                BoxData<double, NCOMP, HOST> boundData(boundBox);
                BoxData<double, NCOMP, HOST> boundErr(boundBox);
                forallInPlace_p(f_phiM, boundData, block, boundX, k, offset);
                patch.copyTo(boundErr);
                boundErr -= boundData;
                boundErr.copyTo(errPatch);
                double e = boundErr.absMax();
                err[nn] = e > err[nn] ? e : err[nn];
            }
        }
#if PR_VERBOSE > 0
        h5.writeMBLevel({"J"}, map, map.jacobian(), "MBInterpOpTests_J_%i", nn);
        h5.writeMBLevel({"phi"}, map, hostSrc, "MBInterpOpTests_Src_%i", nn);
        h5.writeMBLevel({"phi"}, map, hostDst, "MBInterpOpTests_Dst_%i", nn);
        h5.writeMBLevel({"err"}, map, hostErr, "MBInterpOpTests_Err_%i", nn);
        for (auto iter : layout)
        {
            for (auto pi : Box::Kernel(1))
            {
                auto bounds = map.jacobian().bounds(iter, pi);
                pout() << "patch: " << iter.global() << " | dir: " << pi << std::endl;
                for (auto bi : bounds)
                {
                    bi.localData->printData();
                }
            }
        }
#endif
        PR_DEBUG_MSG(1, "Error: %3.2e", err[nn]);
        domainSize *= 2;
        boxSize *= 2;
    }

    double rateTol = 0.1;
    for (int ii = 1; ii < N; ii++)
    {
        double rate = log(err[ii-1]/err[ii])/log(2.0);
        PR_DEBUG_MSG(1,"Rate: %3.2f", rate);
        double rateErr = abs(rate - 4);
        EXPECT_LT(rateErr, rateTol);
    }
}
#endif
int main(int argc, char *argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
#ifdef PR_MPI
    MPI_Init(&argc, &argv);
#endif
    PR_TIMER_SETFILE("MBInterpOpTests.D" + to_string(DIM) + ".N" + to_string(numProc())
        + ".time.table");
    int result = RUN_ALL_TESTS();
    PR_TIMER_REPORT();
#ifdef PR_MPI
    MPI_Finalize();
#endif
    return result;
}
