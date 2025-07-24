#include <gtest/gtest.h>
#include "Proto.H"
#include "TestFunctions.H"
#include "MBMap_Shear.H"
#include "MBMap_XPointRigid.H"

using namespace Proto;

#if 0
TEST(MBInterpOp, MBDataPoint)
{

    // grid parameters
    int domainSize = 32;
    int boxSize = 16;
    int ghostSize = 1;
    
    // interplating function parameters
    double order = 4.0;
    Array<double, DIM> exp{1,1,0,0,0,0};
    exp *= order;
    Array<double, DIM> offset{0,0,0,0,0,0};
    offset += 0.1;
  
    auto domain = buildShear(domainSize);
    Point boxSizeVect = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);
    
    // initialize data and map
    MBLevelBoxData<double, 1, HOST> hostSrc(layout, Point::Ones(ghostSize));
    MBLevelBoxData<double, 1, HOST> hostDst(layout, Point::Ones(ghostSize));

    hostSrc.setVal(1);
    hostDst.setVal(2);

    MBLevelMap<MBMap_Shear<HOST>, HOST> map;
    map.define(layout, Point::Ones(ghostSize));
    map.initialize();
    
    std::vector<MBDataPoint> points;
    for (auto iter : layout)
    {
        for (auto pi : layout[iter])
        {
            points.push_back(MBDataPoint(iter, pi, layout));
        }
    }
#if PR_VERBOSE > 0
    HDF5Handler h5;
    h5.writeMBLevel({"v"},map, hostSrc, "SRC_0"); 
    h5.writeMBLevel({"v"},map, hostDst, "DST_0"); 
#endif
    
    for (auto& pi : points)
    {
        auto vSrc = hostSrc[pi];
        auto vDst = hostDst[pi];
        vDst(0) = vSrc(0);
    }

#if PR_VERBOSE > 0
    h5.writeMBLevel({"v"},map, hostSrc, "SRC_1"); 
    h5.writeMBLevel({"v"},map, hostDst, "DST_1"); 
#endif
}
#endif
#if DIM == 2
#if 0
TEST(MBInterpOp, ShearDebug)
{
#if PR_VERBOSE > 0
    HDF5Handler h5;
#endif
    int domainSize = 16;
    int boxSize = 8;
    int ghostSize = 2;
    
    // initialize data
    auto domain = buildShear(domainSize);
    Point boxSizeVect = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);
    // initialize data and map
    MBLevelMap<MBMap_Shear<HOST>, HOST> map;
    map.define(layout, Point::Ones(ghostSize));
    map.initialize();

    MBInterpOp interp(map);
    Box boundBox = Box::Cube(domainSize).adjacent(Point::X(), ghostSize);
    int ii = 0;
    for (auto iter : layout)
    {
        if (layout.block(iter) != 0) { continue; }
        Box patchBox = layout[iter].grow(ghostSize);
        for (auto pi : patchBox)
        {
            if (boundBox.containsPoint(pi))
            {
                MBDataPoint di(iter, pi, layout);
                auto& pointOp = interp(di);
                string filename = "INTERP_" + std::to_string(ii);
                pointOp.writeFootprint(filename);
                ii++;
            }
        }
    }
    HDF5Handler h5;
    h5.writeMBLevel(map, map.map(), "MAP");
}
#endif
#if 1
TEST(MBInterpOp, ShearTest)
{
#if PR_VERBOSE > 0
    HDF5Handler h5;
#endif
    int domainSize = 16;
    int boxSize = 8;
    int ghostSize = 2;
    
    double order = 4.0;
    Array<double, DIM> exp{1,1,0,0,0,0};
    exp *= order;
    Array<double, DIM> offset{0,0,0,0,0,0};
    offset += 0.1;
  
    int numIter = 3;
    double err[numIter];
    for (int nn = 0; nn < numIter; nn++)
    {
        err[nn] = 0;
        // initialize data
        auto domain = buildShear(domainSize);
        Point boxSizeVect = Point::Ones(boxSize);
        MBDisjointBoxLayout layout(domain, boxSizeVect);
        // initialize data and map
        MBLevelBoxData<double, 1, HOST> hostSrc(layout, Point::Ones(ghostSize));
        MBLevelBoxData<double, 1, HOST> hostDst(layout, Point::Ones(ghostSize));
        MBLevelBoxData<double, 1, HOST> hostErr(layout, Point::Ones(ghostSize));
        MBLevelMap<MBMap_Shear<HOST>, HOST> map;
        map.define(layout, Point::Ones(ghostSize));
        map.initialize();

        auto C2C = Stencil<double>::CornersToCells(4);
        for (auto iter : layout)
        {
            auto block = layout.block(iter);
            auto& src_i = hostSrc[iter];
            auto& dst_i = hostDst[iter];
            Box b_i = C2C.domain(layout[iter]).grow(ghostSize);
            BoxData<double, DIM> x_i(b_i.grow(PR_NODE));
            // Jacobian and NT are computed but not used
            BoxData<double, 1> J_i(b_i);
            FluxBoxData<double, DIM> NT(b_i);
            map.apply(x_i, J_i, NT, block);
            BoxData<double, 1> x_pow = forall_p<double, 1>(f_polyM, block, x_i, exp, offset);
            src_i |= C2C(x_pow);
            dst_i |= C2C(x_pow);
        }
        hostErr.setVal(0);
        hostDst.exchange();
        MBInterpOp interp(map);
#if PR_VERBOSE > 0
        h5.writeMBLevel({"interp"}, map, hostDst, "MBInterpOpTests_Shear_Dst_N%i_0",nn);
#endif
        interp.apply(hostDst, hostDst);
        for (auto iter : layout)
        {
            auto& src_i = hostSrc[iter];
            auto& dst_i = hostDst[iter];
            auto& err_i = hostErr[iter];
            for (auto dir : Box::Kernel(1))
            {
                if (layout.isBlockBoundary(layout.patch(iter), dir))
                {
                    Box boundBox = layout[iter].adjacent(dir*ghostSize);
                    BoxData<double, 1, HOST> error(boundBox);
                    dst_i.copyTo(error);
                    error -= src_i;
                    err[nn] = max(error.absMax(), err[nn]);
                    error.copyTo(err_i);
                }
            }
        }
#if PR_VERBOSE > 0
        std::cout << "Error (Max Norm): " << err[nn] << std::endl;
        h5.writeMBLevel({"soln"}, map, hostSrc, "MBInterpOpTests_Shear_Src_N%i",nn);
        h5.writeMBLevel({"interp"}, map, hostDst, "MBInterpOpTests_Shear_Dst_N%i_1",nn);
        h5.writeMBLevel({"err"}, map, hostErr, "MBInterpOpTests_Shear_Err_N%i",nn);
#endif
        domainSize *= 2;
        boxSize *= 2;
    }
    for (int ii = 1; ii < numIter; ii++)
    {
        double rate = log(err[ii-1]/err[ii])/log(2.0);
#if PR_VERBOSE > 0
        std::cout << "Convergence Rate: " << rate << std::endl;
#endif
        EXPECT_TRUE(err[ii] < 1e-12 || rate > order - 0.5);
    }
}
#endif

namespace {
    template<int NBLOCK>
    double computeXPointInterpError(MBDisjointBoxLayout& layout, int ghostSize, int order, int refIter)
    {
        #if PR_VERBOSE > 0
        HDF5Handler h5;
        #endif
        Array<double, DIM> exp{1.0,1.0,0,0,0,0};
        exp *= order;
        Array<double, DIM> offset{0,0,0,0,0,0};
        offset += 0.1;

        // initialize data and map
        MBLevelBoxData<double, 1, HOST> hostSrc(layout, Point::Ones(ghostSize));
        MBLevelBoxData<double, 1, HOST> hostDst(layout, Point::Ones(ghostSize));
        MBLevelBoxData<double, 1, HOST> hostErr(layout, Point::Ones(ghostSize));

        MBLevelMap<MBMap_XPointRigid<NBLOCK, HOST>, HOST> map;
        map.define(layout, Point::Ones(ghostSize));
        auto C2C = Stencil<double>::CornersToCells(4);
        for (auto iter : layout)
        {
            auto block = layout.block(iter);
            auto &src_i = hostSrc[iter];
            auto &dst_i = hostDst[iter];
            Box b_i = C2C.domain(layout[iter]).grow(ghostSize);
            BoxData<double, DIM> x_i(b_i.grow(PR_NODE));
            // Jacobian and NT are computed but not used
            BoxData<double, 1> J_i(b_i);
            FluxBoxData<double, DIM> NT(b_i);
            map.apply(x_i, J_i, NT, block);
            BoxData<double, 1> x_pow = forall_p<double, 1>(f_polyM, block, x_i, exp, offset);
            src_i |= C2C(x_pow);
            dst_i |= C2C(x_pow);
        }
        hostErr.setVal(0);
#if PR_VERBOSE > 0
        h5.writeMBLevel({"data"}, map, hostDst, "MBInterpOpTests_XPoint_Data_N%i_R%i_0", NBLOCK, refIter);
        h5.writeMBLevelBoundsUnified({"data"}, hostDst, "MBInterpOpTests_XPoint_DataBounds_N%i_R%i_0", NBLOCK, refIter);
#endif
        hostDst.exchange(); // fill boundary data
#if PR_VERBOSE > 0
        h5.writeMBLevelBoundsUnified({"data"}, hostDst, "MBInterpOpTests_XPoint_DataBounds_N%i_R%i_1", NBLOCK, refIter);
#endif
        MBInterpOp interp(map);
        interp.apply(hostDst, hostDst);
#if PR_VERBOSE > 0
        h5.writeMBLevel({"data"}, map, hostDst, "MBInterpOpTests_XPoint_Data_N%i_R%i_1", NBLOCK, refIter);
#endif
        double errNorm = 0;
        for (auto iter : layout)
        {
            auto &src_i = hostSrc[iter];
            auto &dst_i = hostDst[iter];
            auto &err_i = hostErr[iter];
            for (auto dir : Box::Kernel(1))
            {
                if (layout.isBlockBoundary(layout.patch(iter), dir))
                {
                    Box boundBox = layout[iter].adjacent(dir * ghostSize);
                    BoxData<double, 1, HOST> error(boundBox);
                    dst_i.copyTo(error);
                    error -= src_i;
                    errNorm = max(error.absMax(), errNorm);
                    error.copyTo(err_i);
                }
            }
        }

#if PR_VERBOSE > 0
        h5.writeMBLevel({"err"}, map, hostErr, "MBInterpOpTests_XPoint_Err_N%i_%i", NBLOCK, refIter);
#endif
        return errNorm;
    }
    
}
#if 1
TEST(MBInterpOp, XPointTest)
{
#if PR_VERBOSE > 0
    HDF5Handler h5;
#endif
    int domainSize = 16;
    int boxSize = 16;
    int ghostSize = 2;
    int numIter = 3;
    double order = 4;

    double err3[numIter];
    double err4[numIter];
    double err5[numIter];
    for (int nn = 0; nn < numIter; nn++)
    {   
        auto domain3 = buildXPoint(domainSize, 3);
        auto domain4 = buildXPoint(domainSize, 4);
        auto domain5 = buildXPoint(domainSize, 5);
        MBDisjointBoxLayout layout3(domain3, Point::Ones(boxSize));
        MBDisjointBoxLayout layout4(domain4, Point::Ones(boxSize));
        MBDisjointBoxLayout layout5(domain5, Point::Ones(boxSize));
        err3[nn] = computeXPointInterpError<3>(layout3, ghostSize, order, nn);
        err4[nn] = computeXPointInterpError<4>(layout4, ghostSize, order, nn);
        err5[nn] = computeXPointInterpError<5>(layout5, ghostSize, order, nn);
        domainSize *= 2;
    }
    for (int ii = 1; ii < numIter; ii++)
    {
        double rate3 = log(err3[ii - 1] / err3[ii]) / log(2.0);
        double rate4 = log(err4[ii - 1] / err4[ii]) / log(2.0);
        double rate5 = log(err5[ii - 1] / err5[ii]) / log(2.0);
#if PR_VERBOSE > 0
        std::cout << "3 Blocks | Error (Max Norm): " << err3[ii] << " | Convergence Rate: " << rate3 << std::endl;
        std::cout << "4 Blocks | Error (Max Norm): " << err4[ii] << " | Convergence Rate: " << rate4 << std::endl;
        std::cout << "5 Blocks | Error (Max Norm): " << err5[ii] << " | Convergence Rate: " << rate5 << std::endl;
#endif
        EXPECT_TRUE(err3[ii] < 1e-12 || rate3 > order - 0.5);
        EXPECT_TRUE(err4[ii] < 1e-12 || rate4 > order - 0.5);
        EXPECT_TRUE(err5[ii] < 1e-12 || rate5 > order - 0.5);
    }
}
#endif
#if 1
TEST(MBInterpOp, XPointRefined)
{
    int ghostSize = 2;
    int numIter = 2;
    double order = 4;
    int domainSize = 16;
    int boxSize = 8;
    constexpr int numBlocks = 5;
    Array<double, DIM> exp{1.0,1.0,0,0,0,0};
    exp *= order;
    Array<double, DIM> offset{0,0,0,0,0,0};
    offset += 0.1;

    double errNorm[numIter];
    for (int nn = 0; nn < numIter; nn++)
    {
        auto domain = buildXPoint(domainSize, numBlocks);
        std::vector<MBPoint> patches;
        std::vector<Point> boxSizes;
        for (int block = 0; block < numBlocks; block++)
        {
            //TODO: Fix this case (grab more interior data...)
            if (false)
            {
                for (auto pi : Box::Cube(domainSize / boxSize))
                {
                    patches.push_back(MBPoint(pi, block));
                }
            } else {
                patches.push_back(MBPoint(Point::Ones(domainSize / boxSize - 1), block));
            }
            boxSizes.push_back(Point::Ones(boxSize));
        }
        MBDisjointBoxLayout layout(domain, patches, boxSizes);
        // pr_out() << "Test Layout: " << std::endl;
        // layout.print();
        errNorm[nn] = computeXPointInterpError<numBlocks>(layout, ghostSize, 4, nn);
        domainSize *= 2;
    }
    for (int ii = 1; ii < numIter; ii++)
    {
        double rate = log(errNorm[ii-1]/errNorm[ii])/log(2.0);
        #if PR_VERBOSE > 0
        std::cout << "Error (Max Norm): " << errNorm[ii] << " | Convergence Rate: " << rate << std::endl;
        #endif
        EXPECT_NEAR(rate, 4.0, 0.3);
    }
}
#endif

#endif
#if DIM == 3
#if 1
TEST(MBInterpOp, CubedSphereShellTest)
{
#define CUBED_SPHERE_SHELL_R0 0.5
#define CUBED_SPHERE_SHELL_R1 1.0
#if PR_VERBOSE > 0
    HDF5Handler h5;
#endif
    int domainSize = 16;
    int boxSize = 8;
    int thickness = 16;
    int ghostSize = 2;
    bool cullRadialGhost = false;
    double order = 4.0;
    int radialDir = CUBED_SPHERE_SHELL_RADIAL_COORD;
    Array<double, DIM> exp{1,1,1,0,0,0};
    exp *= order;
    Array<double, DIM> offset{0.1,0.2,0.3,0,0,0};
    Point ghost = Point::Ones(ghostSize);
    if (cullRadialGhost) { ghost[radialDir] = 0;}
    int N = 3;
    double err[N];
    double errL1[N];
    for (int nn = 0; nn < N; nn++)
    {
        err[nn] = 0.0;
        errL1[nn] = 0.0;
        auto domain = CubedSphereShell::Domain(domainSize, thickness, radialDir);
        Point boxSizeVect = Point::Ones(boxSize);
        boxSizeVect[radialDir] = min(thickness, boxSize);
        MBDisjointBoxLayout layout(domain, boxSizeVect);

        auto map = CubedSphereShell::Map<HOST>(layout, ghost);
        // initialize data
        MBLevelBoxData<double, 1, HOST> hostSrc(layout, ghost);
        MBLevelBoxData<double, 1, HOST> hostDst(layout, ghost);
        MBLevelBoxData<double, 1, HOST> hostErr(layout, ghost);
        auto C2C = Stencil<double>::CornersToCells(4);
        for (auto iter : layout)
        {
            auto block = layout.block(iter);
            auto& src_i = hostSrc[iter];
            Box b_i = C2C.domain(layout[iter]).grow(ghost);
            BoxData<double, DIM> x_i(b_i.grow(Point::Ones()));
            // Jacobian and NT are computed but not used
            BoxData<double, 1> J_i(x_i.box());
            FluxBoxData<double, DIM> NT(layout[iter]);
            map.apply(x_i, J_i, NT, block);
            BoxData<double, 1> x_pow = forall_p<double, 1>(f_polyM, block, x_i, exp, offset);
            src_i |= C2C(x_pow);
        }

        hostSrc.exchange();
        hostDst.setVal(0);
        hostErr.setVal(0);
    
        MBInterpOp op = CubedSphereShell::InterpOp<HOST>(layout, ghost, order);
        // apply the operator on all block boundaries
        op.apply(hostDst, hostSrc);
        hostDst.exchange();
        for (auto iter : layout)
        {
            auto block = layout.block(iter);
            Box blockDomain = layout.domain().getBlock(block).box();

            auto& err_i = hostErr[iter];
            auto& dst_i = hostDst[iter];
            auto& src_i = hostSrc[iter];

            BoxData<double, 1> tmp(layout[iter]);
            src_i.copyTo(tmp);
            dst_i.copyTo(err_i);
            err_i -= src_i;
            err_i += tmp;
            err_i.setVal(0,blockDomain & err_i.box()); 
            for (auto dir : Box::Kernel(1).grow(-Point::Basis(radialDir)))
            {
                if (dir == Point::Zeros()) {continue; }
                Box bi = blockDomain.adjacent(dir*ghost);
                BoxData<double> ei(bi);
                ei.setVal(0);
                err_i.copyTo(ei);
                err[nn] = max(ei.absMax(), err[nn]);
                errL1[nn] += ei.sum();
            }
        }
        //op.printErrorPoints(hostErr, 1.0);
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
        h5.writeMBLevel({"err"}, map, hostErr, "MBInterpOpTests_CubeSphereShell_Err_%i", nn);
        h5.writeMBLevel({"phi"}, map, hostSrc, "MBInterpOpTests_CubeSphereShell_Src_%i", nn);
        h5.writeMBLevel({"phi"}, map, hostDst, "MBInterpOpTests_CubeSphereShell_Dst_%i", nn);
#endif
        domainSize *= 2;
        boxSize *= 2;
        thickness *= 2;
    }
    for (int ii = 1; ii < N; ii++)
    {
        double rate = log(err[ii-1]/err[ii])/log(2.0);
        double rateL1 = log(errL1[ii-1]/errL1[ii])/log(2.0);
        //EXPECT_GT(rate, 3.5);
        EXPECT_TRUE(err[ii] < 1e-12 || rate >  order - 1.0);
#if PR_VERBOSE > 0
        if (procID() == 0)
        {
            std::cout << "Convergence Rate (Max Norm): " << rate;
            std::cout << " | (L1 Norm): " << rateL1 << std::endl;
        }
        std::cout << "Convergence pass / fail is determined by Max Norm" << std::endl;
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
    PR_TIMER_SETFILE("MBInterpOpTests_DIM" + to_string(DIM) + "_NProc" + to_string(numProc())
            + ".time.table");
    int result = RUN_ALL_TESTS();
    PR_TIMER_REPORT();
#ifdef PR_MPI
    MPI_Finalize();
#endif
    return result;
}
