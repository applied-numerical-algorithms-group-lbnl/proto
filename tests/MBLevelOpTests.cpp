#include <gtest/gtest.h>
#include "Proto.H"
#include "Lambdas.H"

#include "MBLevelMap_CubeSphereShell.H"
#include "MBLevelMap_CubeSphereShellPolar.H"
#include "MBLevelMap_XPointRigid.H"
#include "MBLevelMap_Shear.H"
#include "BoxOp_MBLaplace.H"

using namespace Proto;

#if DIM==2
#if 1
TEST(MBLevelOp, ShearLaplace) {

    HDF5Handler h5;
    int domainSize = 32;
    int boxSize = 32;
    int numGhost = 4;
    Array<double, DIM> k{1,1,1,0,0,0};
    Array<double, DIM> offset{0,0,0,0,0,0};
    offset += 0.1;
    Array<Point, DIM+1> srcGhost;
    Array<Point, DIM+1> dstGhost;
    Array<Point, DIM+1> mapGhost;

    srcGhost.fill(Point::Ones(numGhost));
    mapGhost.fill(Point::Ones(numGhost+1));
    dstGhost.fill(Point::Zeros());

    int N = 3;
    double err[N];
    for (int nn = 0; nn < N; nn++)
    {
        err[nn] = 0.0;
        auto domain = buildShear(domainSize);
        Point boxSizeVect = Point::Ones(boxSize);
        MBDisjointBoxLayout layout(domain, boxSizeVect);

        MBLevelMap_Shear<HOST> map;
        map.define(layout, mapGhost);
        
        MBLevelBoxData<double, 1, HOST> hostSrc(layout, srcGhost);
        MBLevelBoxData<double, 1, HOST> hostDst(layout, dstGhost);
        MBLevelBoxData<double, 1, HOST> hostSln(layout, dstGhost);
        MBLevelBoxData<double, 1, HOST> hostErr(layout, dstGhost);

        auto C2C = Stencil<double>::CornersToCells(4);
        for (auto iter : layout)
        {
            auto block = layout.block(iter);
            auto& src_i = hostSrc[iter];
            Box b_i = C2C.domain(layout[iter]).grow(numGhost);
            BoxData<double, DIM> x_i(b_i);
            BoxData<double, 1> J_i(b_i);
            map.apply(x_i, J_i, block);
            BoxData<double, 1> phi = forall_p<double, 1>(f_phiM, block, x_i, k, offset);
            src_i |= C2C(phi);
       
            BoxData<double, 1> lphi = forall_p<double, 1>(f_LphiM, block, x_i, k, offset);
            auto& sln_i = hostSln[iter];
            sln_i |= C2C(lphi);
        }

#if PR_VERBOSE > 0
        h5.writeMBLevel({"phi"}, map, hostSrc, "Shear_Phi_%i", nn);
        h5.writeMBLevel({"Lphi"}, map, hostSln, "Shear_Sln_%i", nn);
#endif

        hostSrc.exchange();
        hostSrc.fillBoundaries();
        hostDst.setVal(0);
        hostErr.setVal(0);
       
        MBLevelOp<BoxOp_MBLaplace, double> op;
        op.define(map);
        op(hostDst, hostSrc);
        hostDst.exchange();
        hostDst.fillBoundaries();
        for (auto iter : layout)
        {
            auto& err_i = hostErr[iter];
            auto& dst_i = hostDst[iter];
            auto& sln_i = hostSln[iter];
            dst_i.copyTo(err_i);
            err_i -= sln_i;
            err[nn] = max(err_i.absMax(), err[nn]);
        }
        Reduction<double, Max> rxn;
        rxn.reduce(&err[nn], 1);
        err[nn] = rxn.fetch();

#if PR_VERBOSE > 0
        if (procID() == 0)
        {
            std::cout << "Error (Max Norm): " << err[nn] << std::endl;
        }
        h5.writeMBLevel({"err"}, map, hostErr, "Shear_Err_%i", nn);
        h5.writeMBLevel({"Lphi"}, map, hostDst, "Shear_LPhi_%i", nn);
        h5.writeMBLevel({"J"}, map, map.jacobian(), "Shear_J_%i", nn);
#endif
        domainSize *= 2;
        boxSize *= 2;
    }
    for (int ii = 1; ii < N; ii++)
    {
        double rate = log(err[ii-1]/err[ii])/log(2.0);
        EXPECT_GT(rate, 3.5);
#if PR_VERBOSE > 0
        if (procID() == 0)
        {
            std::cout << "Convergence Rate (Max Norm): " << rate << std::endl;
        }
#endif
    }
}
#endif
#endif

#if 1
TEST(MBLevelOp, XPointLaplace) {

    HDF5Handler h5;
    int domainSize = 32;
    int boxSize = 32;
    int numGhost = 4;
    int numBlocks = 8;
    Array<double, DIM> k{1,1,1,0,0,0};
    Array<double, DIM> offset{0,0,0,0,0,0};
    offset += 0.1;
    Array<Point, DIM+1> srcGhost;
    Array<Point, DIM+1> dstGhost;
    Array<Point, DIM+1> mapGhost;

    srcGhost.fill(Point::Ones(numGhost));
    mapGhost.fill(Point::Ones(numGhost+1));
    dstGhost.fill(Point::Zeros());

    int N = 3;
    double err[N];
    for (int nn = 0; nn < N; nn++)
    {
        err[nn] = 0.0;
        auto domain = buildXPoint(domainSize, numBlocks);
        Point boxSizeVect = Point::Ones(boxSize);
        MBDisjointBoxLayout layout(domain, boxSizeVect);

        MBLevelMap_XPointRigid<HOST> map;
        map.setNumBlocks(numBlocks);
        map.define(layout, mapGhost);
        
        MBLevelBoxData<double, 1, HOST> hostSrc(layout, srcGhost);
        MBLevelBoxData<double, 1, HOST> hostDst(layout, dstGhost);
        MBLevelBoxData<double, 1, HOST> hostSln(layout, dstGhost);
        MBLevelBoxData<double, 1, HOST> hostErr(layout, dstGhost);

        auto C2C = Stencil<double>::CornersToCells(4);
        for (auto iter : layout)
        {
            auto block = layout.block(iter);
            auto& src_i = hostSrc[iter];
            Box b_i = C2C.domain(layout[iter]).grow(numGhost);
            BoxData<double, DIM> x_i(b_i);
            BoxData<double, 1> J_i(b_i);
            map.apply(x_i, J_i, block);
            BoxData<double, 1> phi = forall_p<double, 1>(f_phiM, block, x_i, k, offset);
            src_i |= C2C(phi);
       
            BoxData<double, 1> lphi = forall_p<double, 1>(f_LphiM, block, x_i, k, offset);
            auto& sln_i = hostSln[iter];
            sln_i |= C2C(lphi);
        }

#if PR_VERBOSE > 0
        h5.writeMBLevel({"phi"}, map, hostSrc, "XPoint_Phi_%i", nn);
        h5.writeMBLevel({"Lphi"}, map, hostSln, "XPoint_Sln_%i", nn);
#endif

        hostSrc.exchange();
        hostSrc.fillBoundaries();
        hostDst.setVal(0);
        hostErr.setVal(0);
       
        MBLevelOp<BoxOp_MBLaplace, double> op;
        op.define(map);
        op(hostDst, hostSrc);
        hostDst.exchange();
        hostDst.fillBoundaries();
        for (auto iter : layout)
        {
            auto& err_i = hostErr[iter];
            auto& dst_i = hostDst[iter];
            auto& sln_i = hostSln[iter];
            double J0 = map.jacobian()[iter].absMax(); //J is a constant
            dst_i /= (J0);
            dst_i.copyTo(err_i);
            err_i -= sln_i;
            err[nn] = max(err_i.absMax(), err[nn]);
        }
        Reduction<double, Max> rxn;
        rxn.reduce(&err[nn], 1);
        err[nn] = rxn.fetch();

#if PR_VERBOSE > 0
        if (procID() == 0)
        {
            std::cout << "Error (Max Norm): " << err[nn] << std::endl;
        }
        h5.writeMBLevel({"err"}, map, hostErr, "XPoint_Err_%i", nn);
        h5.writeMBLevel({"Lphi"}, map, hostDst, "XPoint_LPhi_%i", nn);
        h5.writeMBLevel({"J"}, map, map.jacobian(), "XPoint_J_%i", nn);
#endif
        domainSize *= 2;
        boxSize *= 2;
    }
    for (int ii = 1; ii < N; ii++)
    {
        double rate = log(err[ii-1]/err[ii])/log(2.0);
        EXPECT_GT(rate, 3.5);
#if PR_VERBOSE > 0
        if (procID() == 0)
        {
            std::cout << "Convergence Rate (Max Norm): " << rate << std::endl;
        }
#endif
    }
}
#endif

#if DIM==3
#if 0
TEST(MBLevelOp, CubeSphereLaplace) {

    HDF5Handler h5;
    int domainSize = 16;
    int boxSize = 16;
    int thickness = 1;
    bool cullRadialGhost = false;
    bool use2DFootprint = true;
    int numGhost = 5;
    int radialDir = CUBE_SPHERE_SHELL_RADIAL_COORD;
    Array<double, DIM> k{1,1,1,0,0,0};
    Array<double, DIM> offset{0,0,0,0,0,0};
    //offset += 0.1;
    Array<Point, DIM+1> dataGhost;
    Array<Point, DIM+1> errGhost;

    dataGhost.fill(Point::Ones(numGhost+2));
    dataGhost[0] = Point::Ones(numGhost);
    if (cullRadialGhost) { dataGhost[0][radialDir] = 0;}
    errGhost.fill(Point::Zeros());

    std::vector<Point> footprint;
    for (auto pi : Box::Kernel(3))
    {
        if (pi.abs().sum() <= 2)
        {
            if (use2DFootprint && (pi[radialDir] != 0)) { continue; }
            footprint.push_back(pi);
        }
    }
    int N = 1;
    double err[N];
    double errL1[N];
    for (int nn = 0; nn < N; nn++)
    {
        err[nn] = 0.0;
        errL1[nn] = 0.0;
        auto domain = buildCubeSphereShell(domainSize, thickness, radialDir);
        Point boxSizeVect = Point::Ones(boxSize);
        boxSizeVect[radialDir] = thickness;
        MBDisjointBoxLayout layout(domain, boxSizeVect);

        MBLevelMap_CubeSphereShell<HOST> map;
        map.define(layout, dataGhost);
        MBLevelMap_CubeSphereShellPolar<HOST> polarMaps[6];
        for (int bi = 0; bi < 6; bi++)
        {
            polarMaps[bi].define(layout, dataGhost, bi);
        }

        MBLevelBoxData<double, 1, HOST> hostSrc(layout, dataGhost);
        MBLevelBoxData<double, 1, HOST> hostDst(layout, dataGhost);
        MBLevelBoxData<double, 1, HOST> hostSln(layout, errGhost);
        MBLevelBoxData<double, 1, HOST> hostErr(layout, errGhost);
        auto C2C = Stencil<double>::CornersToCells(4);
        for (auto iter : layout)
        {
            auto block = layout.block(iter);
            auto& src_i = hostSrc[iter];
            Box b_i = C2C.domain(layout[iter]).grow(dataGhost[0]);
            BoxData<double, DIM> x_i(b_i.grow(Point::Ones()));
            BoxData<double, 1> J_i(layout[iter].grow(Point::Ones() + dataGhost[0]));
            FluxBoxData<double, DIM> NT(layout[iter]);
            map.apply(x_i, J_i, NT, block);
            BoxData<double, 1> phi = forall_p<double, 1>(f_phiM, block, x_i, k, offset);
            src_i |= C2C(phi);
        
            auto& sln_i = hostSln[iter];
            BoxData<double, 1> lphi = forall_p<double, 1>(f_LphiM, block, x_i, k, offset);
            sln_i |= C2C(lphi);
        }

        hostSrc.exchange();
        hostSrc.fillBoundaries();
        hostDst.setVal(0);
        hostErr.setVal(0);
       
        MBLevelOp<BoxOp_MBLaplace, double> op;
        op.define(map);
        op(hostDst, hostSrc);
        for (auto iter : layout)
        {
            auto& err_i = hostErr[iter];
            auto& dst_i = hostDst[iter];
            auto& sln_i = hostSln[iter];
            const auto& J_i = map.jacobian()[iter];
            Operator::cellQuotient(err_i,dst_i, J_i);
            dst_i.copyTo(err_i);
            err_i -= sln_i;
            err[nn] = max(err_i.absMax(), err[nn]);
            errL1[nn] += err_i.sum();
        }
        Reduction<double, Max> rxn;
        rxn.reduce(&err[nn], 1);
        err[nn] = rxn.fetch();
        rxn.reset();
        rxn.reduce(&errL1[nn], 1);
        errL1[nn] = rxn.fetch() / domainSize;

#if PR_VERBOSE > 0
        if (procID() == 0)
        {
            std::cout << "Error (Max Norm): " << err[nn];
            std::cout << " | Error (L1 Norm): " << errL1[nn] << std::endl;
        }
        h5.writeMBLevel({"err"}, map, hostErr, "CubeSphereShell_Err_%i", nn);
        h5.writeMBLevel({"phi"}, map, hostSrc, "CubeSphereShell_Phi_%i", nn);
        h5.writeMBLevel({"Lphi"}, map, hostDst, "CubeSphereShell_LPhi_%i", nn);
        h5.writeMBLevel({"Lphi"}, map, hostSln, "CubeSphereShell_Sln_%i", nn);
        h5.writeMBLevel({"J"}, map, map.jacobian(), "CubeSphereShell_J%i", nn);
#endif
        domainSize *= 2;
        boxSize *= 2;
    }
    for (int ii = 1; ii < N; ii++)
    {
        double rate = log(err[ii-1]/err[ii])/log(2.0);
        double rateL1 = log(errL1[ii-1]/errL1[ii])/log(2.0);
        EXPECT_GT(rate, 3.5);
        EXPECT_GT(rateL1, 3.5);
#if PR_VERBOSE > 0
        if (procID() == 0)
        {
            std::cout << "Convergence Rate (Max Norm): " << rate;
            std::cout << " | (L1 Norm): " << rateL1 << std::endl;
        }
#endif
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
