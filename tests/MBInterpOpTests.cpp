#include <gtest/gtest.h>
#include "Proto.H"
#include "Lambdas.H"
#include "MBLevelMap_Shear.H"
#include "MBLevelMap_CubeSphereShell.H"

using namespace Proto;
#if DIM == 2
TEST(MBInterpOp, ShearTest)
{
    int domainSize = 4;
    int boxSize = 4;
    Array<double, DIM> exp{4,4,0,0,0,0};
    Array<double, DIM> offset{0,0,0,0,0,0};
    HDF5Handler h5;
    Array<Point, DIM+1> ghost;
    ghost.fill(Point::Ones(4));
    ghost[0] = Point::Ones(1);
    std::vector<Point> footprint;
    for (auto pi : Box::Kernel(3))
    {
        if (pi.abs().sum() <= 2)
        {
            footprint.push_back(pi);
        }
    }
    int N = 3;
    double err[N];
    for (int nn = 0; nn < N; nn++)
    {
        err[nn] = 0.0;
        auto domain = buildShear(domainSize);
        Point boxSizeVect = Point::Ones(boxSize);
        MBDisjointBoxLayout layout(domain, boxSizeVect);

        ghost[0] = Point::Ones(2);
        MBLevelMap_Shear<HOST> map;
        map.define(layout, ghost);

        ghost[0] = Point::Ones(1);
        MBLevelBoxData<double, 1, HOST> hostSrc(layout, ghost);
        MBLevelBoxData<double, 1, HOST> hostDst(layout, ghost);
        MBLevelBoxData<double, 1, HOST> hostErr(layout, ghost);

        for (auto iter : layout)
        {
            auto& src_i = hostSrc[iter];
            auto& x_i = map.map()[iter];
            auto block = layout.block(iter);
            BoxData<double, 1> x_pow = forall_p<double, 1>(f_polyM, block, x_i, exp, offset);
            src_i |= Stencil<double>::CornersToCells(4)(x_pow);
        }

        //hostSrc.initialize(f_polyM, map.map(), dx, exp, offset);
        hostSrc.exchange();
        hostSrc.fillBoundaries();
        hostDst.setVal(0);
        hostErr.setVal(0);
        

        auto blockDomainBox = Box::Cube(domainSize);
        for (auto iter : layout)
        {
            auto block = layout.block(iter);
            Box patchBox = layout[iter];
            for (auto dir : Box::Kernel(1))
            {
                auto bounds = hostSrc.bounds(iter, dir);
                for (auto bound : bounds)
                {
                    Box boundBox = patchBox.adjacent(ghost[0]*dir);
                    if (blockDomainBox.contains(boundBox)) { continue; }
                    for (auto bi : boundBox)
                    {
                        MBDataPoint dstDataPoint(iter, bi, layout);
                        MBPointInterpOp op(dstDataPoint, ghost[0], map, footprint, 4);
                        op.apply(hostDst, hostSrc);
#if PR_VERBOSE > 1
                        pout() << "Coefs at point " << bi << std::endl;
                        auto coefs = op.coefs(hostSrc);
                        coefs.print("%10.2e");
#endif
                        double interpValue = hostDst[dstDataPoint](0);
                        double exactValue =  hostSrc[dstDataPoint](0);
                        double errorValue = interpValue - exactValue;
                        err[nn] = max(abs(errorValue), err[nn]);
                        hostErr[dstDataPoint](0) = std::abs(errorValue);
                    }
                }
            }
        }
#if PR_VERBOSE > 0
        std::cout << "Error (Max Norm): " << err[nn] << std::endl;
        h5.writeMBLevel({"err"}, map, hostErr, "MBInterpOpTests_Shear_Err_%i", nn);
        h5.writeMBLevel({"phi"}, map, hostSrc, "MBInterpOpTests_Shear_Src_%i", nn);
        h5.writeMBLevel({"phi"}, map, hostDst, "MBInterpOpTests_Shear_Dst_%i", nn);
#endif
        domainSize *= 2;
        boxSize *= 2;
    }

    for (int ii = 1; ii < N; ii++)
    {
        double rate = log(err[ii-1]/err[ii])/log(2.0);
#if PR_VERBOSE > 0
        std::cout << "Convergence Rate: " << rate << std::endl;
#endif
    }
}
#endif
#if DIM == 3
TEST(MBInterpOp, CubeSphereShellTest)
{
    HDF5Handler h5;
    int domainSize = 4;
    int boxSize = 4;
    int thickness = 1;
    int radialDir = CUBE_SPHERE_SHELL_RADIAL_COORD;
    Array<double, DIM> exp{0,0,0,0,0,0};
    Array<double, DIM> offset{0,0,0,0,0,0};
    Array<Point, DIM+1> ghost;
    ghost.fill(Point::Ones(4));
    ghost[0] = Point::Ones(1);
    ghost[0][radialDir] = 0;
    std::vector<Point> footprint;
    for (auto pi : Box::Kernel(3))
    {
        if (pi.abs().sum() <= 2 && (pi[radialDir] == 0))
        {
            footprint.push_back(pi);
        }
    }
    int N = 3;
    double err[N];
    for (int nn = 0; nn < N; nn++)
    {
        err[nn] = 0.0;
        auto domain = buildCubeSphereShell(domainSize, thickness, radialDir);
        Point boxSizeVect = Point::Ones(boxSize);
        boxSizeVect[radialDir] = thickness;
        MBDisjointBoxLayout layout(domain, boxSizeVect);

        ghost[0] = Point::Ones(3);
        ghost[0][radialDir] = 0;
        MBLevelMap_CubeSphereShell<HOST> map;
        map.define(layout, ghost);

        ghost[0] = Point::Ones(1);
        ghost[0][radialDir] = 0;

        MBLevelBoxData<double, 1, HOST> hostSrc(layout, ghost);
        MBLevelBoxData<double, 1, HOST> hostDst(layout, ghost);
        MBLevelBoxData<double, 1, HOST> hostErr(layout, ghost);

        for (auto iter : layout)
        {
            auto& src_i = hostSrc[iter];
            auto& x_i = map.map()[iter];
            auto block = layout.block(iter);
            BoxData<double, 1> x_pow = forall_p<double, 1>(f_polyM, block, x_i, exp, offset);
            src_i |= Stencil<double>::CornersToCells(4)(x_pow);
        }

        //hostSrc.initialize(f_polyM, map.map(), dx, exp, offset);
        hostSrc.exchange();
        hostSrc.fillBoundaries();
        hostDst.setVal(0);
        hostErr.setVal(0);
        

        auto blockDomainBox = Box::Cube(domainSize);
        for (auto iter : layout)
        {
            auto block = layout.block(iter);
            Box patchBox = layout[iter];
            for (auto dir : Box::Kernel(1))
            {
                auto bounds = hostSrc.bounds(iter, dir);
                for (auto bound : bounds)
                {
                    Box boundBox = patchBox.adjacent(ghost[0]*dir);
                    if (blockDomainBox.contains(boundBox)) { continue; }
                    for (auto bi : boundBox)
                    {
                        MBDataPoint dstDataPoint(iter, bi, layout);
                        MBPointInterpOp op(dstDataPoint, ghost[0], map, footprint, 4);
                        op.apply(hostDst, hostSrc);
#if PR_VERBOSE > 1
                        pout() << "Coefs at point " << bi << std::endl;
                        auto coefs = op.coefs(hostSrc);
                        coefs.print("%10.2e");
#endif
                        double interpValue = hostDst[dstDataPoint](0);
                        double exactValue =  hostSrc[dstDataPoint](0);
                        double errorValue = interpValue - exactValue;
                        err[nn] = max(abs(errorValue), err[nn]);
                        hostErr[dstDataPoint](0) = std::abs(errorValue);
                    }
                }
            }
        }
#if PR_VERBOSE > 0
        std::cout << "Error (Max Norm): " << err[nn] << std::endl;
        h5.writeMBLevel({"err"}, map, hostErr, "MBInterpOpTests_CubeSphereShell_Err_%i", nn);
        h5.writeMBLevel({"phi"}, map, hostSrc, "MBInterpOpTests_CubeSphereShell_Src_%i", nn);
        h5.writeMBLevel({"phi"}, map, hostDst, "MBInterpOpTests_CubeSphereShell_Dst_%i", nn);
#endif
        domainSize *= 2;
        boxSize *= 2;
    }

    for (int ii = 1; ii < N; ii++)
    {
        double rate = log(err[ii-1]/err[ii])/log(2.0);
#if PR_VERBOSE > 0
        std::cout << "Convergence Rate: " << rate << std::endl;
#endif
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
