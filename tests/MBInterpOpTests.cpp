#include <gtest/gtest.h>
#include "Proto.H"
#include "Lambdas.H"

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
void computeC(MBMap<MAP>& a_map, int a_boxSize, Point a_dst, Point a_src)
{
    double h = 1.0/a_boxSize;
    auto C2C = Stencil<double>::cornersToCells(4);
    
    Box Bx = Box::Cube(a_boxSize*2);
    Box Bx0 = C2C.domain(Bx).extrude(Point::Ones(),-1);
   
    std::vector<BoxData<double, DIM>> X;
    auto Xg = a_map(a_dst, Point::Zeros(), 0);
    for (int ii = 0; ii < 4; ii++)
    {
        X.push_back(a_map(Bx0,ii,0));
        //X[ii] -= Xg;
    }

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
  
            pout() << "X ^ " << exponents[ei] << " | block: " << bi << std::endl;
            Xp[bi][ei].printData(4);
            Xp_avg[bi][ei].printData(4);
        }
    }
}

TEST(MBPointInterpOp, CheckMatrix)
{

    auto C2C = Stencil<double>::cornersToCells(4);

    int domainSize = 4;
    int boxSize = 4;
    
    auto domain = buildShear(domainSize);
    Point boxSizeVect = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);

    Array<Point, DIM+1> ghost;
    ghost.fill(Point::Ones(4));
    ghost[0] = Point::Ones(1);
    Point boundGhost = Point::Ones();

    // initialize map
    MBMap<ShearMap_t> map(ShearMap, layout, ghost, boundGhost);

}

#if 0
TEST(MBPointInterpOp, ShearApply) {
    int domainSize = 8;
    int boxSize = 8;
    HDF5Handler h5;

    constexpr int N = 3;
    Array<double, N> err(0);
    for (int nn = 0; nn < N; nn++)
    {
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

        Array<double, DIM> exp{0,2,0,0,0,0};
        Array<double, DIM> offset{0,0,0,0,0,0};
        hostSrc.initialize(f_polyM, map, exp, offset);
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
                auto boundX = map(boundBox, block, PR_CELL);
                BoxData<double, NCOMP, HOST> boundData(boundBox);
                BoxData<double, NCOMP, HOST> boundErr(boundBox);
                forallInPlace_p(f_polyM, boundData, block, boundX, exp, offset);
                patch.copyTo(boundErr);
                boundErr -= boundData;
                boundErr.copyTo(errPatch);
                double e = boundErr.absMax();
                err[nn] = e > err[nn] ? e : err[nn];
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
