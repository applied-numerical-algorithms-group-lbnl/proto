#include <gtest/gtest.h>
#include "Proto.H"
#include "Lambdas.H"

#include "MBMap_XPointRigid.H"
#include "MBMap_Shear.H"
#include "BoxOp_MBLaplace.H"

using namespace Proto;



#if DIM==2
#if 1
TEST(MBLevelOp, Iteration) {
    HDF5Handler h5;
    int domainSize = 16;
    int boxSize = 16;
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

    auto domain = buildShear(domainSize);
    Point boxSizeVect = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);
    
    MBLevelMap<MBMap_Shear, HOST> map;
    map.define(layout, mapGhost);

    MBLevelOp<BoxOp_MBLaplace, MBMap_Shear, double> op(map);

    for (auto iter : layout)
    {
        EXPECT_EQ(op[iter].box(), layout[iter]);
    }
}
#endif
#if 1
TEST(MBLevelOp, ShearLaplace) {

    HDF5Handler h5;
    int domainSize = 16;
    int boxSize = 16;
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

        MBLevelBoxData<double, 1, HOST> hostSrc(layout, srcGhost);
        MBLevelBoxData<double, 1, HOST> hostDst(layout, dstGhost);
        MBLevelBoxData<double, 1, HOST> hostSln(layout, dstGhost);
        MBLevelBoxData<double, 1, HOST> hostErr(layout, dstGhost);
        
        MBLevelMap<MBMap_Shear, HOST> map;
        map.define(layout, mapGhost);

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
        hostDst.setVal(0);
        hostErr.setVal(0);
       
        MBLevelOp<BoxOp_MBLaplace, MBMap_Shear, double> op;
        op.define(map);
        op(hostDst, hostSrc);
        hostDst.exchange();
        for (auto iter : layout)
        {
            auto& err_i = hostErr[iter];
            auto& dst_i = hostDst[iter];
            auto& sln_i = hostSln[iter];
            dst_i.copyTo(err_i);
            double J0 = map.jacobian()[iter].absMax(); //J is a constant
            err_i /= J0;
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
#if DIM==2
TEST(MBLevelOp, XPointLaplace) {

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
        double dx = 1.0/domainSize;
        auto domain = buildXPoint(domainSize);
        Point boxSizeVect = Point::Ones(boxSize);
        MBDisjointBoxLayout layout(domain, boxSizeVect);

        MBLevelMap<MBMap_XPointRigid, HOST> map;
        map.define(layout, mapGhost);
        
        MBLevelBoxData<double, 1, HOST> hostSrc(layout, srcGhost);
        MBLevelBoxData<double, DIM, HOST> hostFlx(layout, srcGhost);
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
            double J0 = J_i.absMax(); //J is a constant
            BoxData<double, 1> phi = forall_p<double, 1>(f_phiM, block, x_i, k, offset);
            src_i |= C2C(phi);
       
            BoxData<double, 1> lphi = forall_p<double, 1>(f_LphiM, block, x_i, k, offset);
            auto& sln_i = hostSln[iter];
            sln_i |= C2C(lphi);
            sln_i *= J0;
        }

#if PR_VERBOSE > 0
        h5.writeMBLevel({"phi"}, map, hostSrc, "XPoint_Phi_%i", nn);
        h5.writeMBLevel({"Lphi"}, map, hostSln, "XPoint_Sln_%i", nn);
#endif

        hostSrc.exchange();
        hostDst.setVal(0);
        hostErr.setVal(0);
        hostFlx.setVal(7);

        MBLevelOp<BoxOp_MBLaplace, MBMap_XPointRigid, double> op;
        op.define(map);
        op(hostDst, hostSrc);
        hostDst.exchange();
         
        for (auto iter : layout)
        {
            auto& src_i = hostSrc[iter];    //source data already initialized
            auto& flx_i = hostFlx[iter];    //uninitialized data with DIM components
            for (int ii = 0; ii < DIM; ii++)
            {
                auto fd = slice(flx_i, ii); //alias to a single component(?)
                op[iter].flux(fd, src_i,ii);              //update fd
            }

            auto& err_i = hostErr[iter];
            auto& dst_i = hostDst[iter];
            auto& sln_i = hostSln[iter];
            //double J0 = map.jacobian()[iter].absMax(); //J is a constant
            //dst_i /= (J0);
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
        h5.writeMBLevel({"F"}, map, hostFlx, "XPoint_Flux_%i", nn);
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

#if 0
#if DIM==2
TEST(MBLevelOp, FluxMatching) {

    HDF5Handler h5;
    int domainSize = 16;
    int boxSize = 16;
    int numGhost = 4;

    auto domain = buildXPoint(domainSize);
    Point boxSizeVect = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);

    MBLevelMap<MBMap_XPointRigid, HOST> map;
    map.define(layout, Point::Ones(numGhost));
    
    MBLevelBoxData<double, 1, HOST> hostSrc(layout, Point::Ones(numGhost));
    MBLevelBoxData<double, 1, HOST> hostDst(layout, Point::Zeros());

    auto C2C = Stencil<double>::CornersToCells(4);
    for (auto iter : layout)
    {
        auto block = layout.block(iter);
        auto& src_i = hostSrc[iter];
        Box b_i = C2C.domain(layout[iter]).grow(numGhost);
        BoxData<double, DIM> x_i(b_i);
        BoxData<double, 1> J_i(b_i);
        map.apply(x_i, J_i, block);
        double J0 = J_i.absMax(); //J is a constant

        BoxData<double, 1> phi = forall_p<double, 1>(
            [] PROTO_LAMBDA (
                Point& a_pt,
                Var<double, 1, HOST>& a_data,
                int a_block,
                Var<double, DIM, HOST>& a_X)
            {
                int coord = a_block % DIM;
                a_data(0) = pow(a_X(coord), 1);
            }, block, x_i);

        src_i |= C2C(phi);
    }
    hostDst.setVal(0);
#if PR_VERBOSE > 0
    h5.writeMBLevel({"phi"}, map, hostSrc, "FluxMatch_Phi_%i",0);
#endif
    MBLevelOp<BoxOp_MBLaplace, MBMap_XPointRigid, double> op;
    op.define(map);
    op(hostDst, hostSrc);
#if PR_VERBOSE > 0
    h5.writeMBLevel({"phi"}, map, hostSrc, "FluxMatch_Phi_%i",1);
    h5.writeMBLevel({"lphi"}, map, hostDst, "FluxMatch_LPhi_%i",0);
#endif
    op.matchFlux(hostDst, hostSrc);
#if PR_VERBOSE > 0
    h5.writeMBLevel({"phi"}, map, hostSrc, "FluxMatch_Phi_%i",2);
    h5.writeMBLevel({"lphi"}, map, hostDst, "FluxMatch_LPhi_%i",1);
#endif
}
#endif
#endif

#if DIM==3
#if 1
TEST(MBLevelOp, CubeSphereLaplace) {

    HDF5Handler h5;
    int domainSize = 16;
    int boxSize = 16;
    int thickness = 1;
    bool cullRadialGhost = false;
    bool use2DFootprint = true;
    int numGhost = 5;
    int radialDir = CUBED_SPHERE_SHELL_RADIAL_COORD;
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
    int N = 2;
    double err[N];
    double errL1[N];
    for (int nn = 0; nn < N; nn++)
    {
        err[nn] = 0.0;
        errL1[nn] = 0.0;
        auto domain = CubedSphereShell::Domain(domainSize, thickness, radialDir);
        Point boxSizeVect = Point::Ones(boxSize);
        boxSizeVect[radialDir] = thickness;
        MBDisjointBoxLayout layout(domain, boxSizeVect);

        MBLevelBoxData<double, 1, HOST> hostSrc(layout, dataGhost);
        MBLevelBoxData<double, 1, HOST> hostDst(layout, dataGhost);
        MBLevelBoxData<double, 1, HOST> hostSln(layout, errGhost);
        MBLevelBoxData<double, 1, HOST> hostErr(layout, errGhost);
        
        auto map = CubedSphereShell::Map(hostSrc); 

        hostSrc.setRandom(0,1);
        h5.writeMBLevel({"phi"}, map, hostSrc, "CubeSphereShell_Phi_Random_%i", nn);
        
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
        hostDst.setVal(0);
        hostErr.setVal(0);
       
        auto op = CubedSphereShell::Operator<BoxOp_MBLaplace, double, HOST>(map);
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
