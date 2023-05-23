#include <gtest/gtest.h>
#include "Proto.H"
#include "Lambdas.H"
#include <iomanip>

#define NCOMP 1

using namespace Proto;
#if 0
TEST(MBPointInterpOp, Shear) {
#if PR_VERBOSE > 0
    pout() << "SHEAR TEST START" << std::endl;
#endif
    int domainSize = 8;
    int boxSize = 8;
    HDF5Handler h5;

    Array<double, DIM> exp{4,4,4,0,0,0};
    Array<double, DIM> offset{0,0,0,0,0,0};
    Array<double, DIM> k{1,1,1,1,1,1};

    // input footprint
    std::vector<Point> footprint;
    for (auto pi : Box::Kernel(3))
    {
        if (pi.abs().sum() <= 2)
        {
            footprint.push_back(pi);
        }
    }

    constexpr int N = 3;
    double err[N];
    for (int nn = 0; nn < N; nn++)
    {
#if PR_VERBOSE > 0
        pout() << "BEGIN_REFINEMENT: domainSize = " << domainSize << std::endl;
#endif
        auto domain = buildShear(domainSize);
        Point boxSizeVect = Point::Ones(boxSize);
        MBDisjointBoxLayout layout(domain, boxSizeVect);

        Array<Point, DIM+1> ghost;
        ghost.fill(Point::Ones(4));
        ghost[0] = Point::Ones(1);
        Point boundGhost = Point::Ones(1);

        // initialize map
        MBMap<ShearMap_t> map(ShearMap, layout, ghost, boundGhost);

        // initialize data
        MBLevelBoxData<double, NCOMP, HOST> hostSrc(layout, ghost);
        MBLevelBoxData<double, NCOMP, HOST> hostDst(layout, ghost);
        MBLevelBoxData<double, NCOMP, HOST> hostErr(layout, ghost);

        hostSrc.initConvolve(f_polyM, map, exp, offset);

        hostSrc.fillBoundaries();
        hostDst.setVal(0);
        hostErr.setVal(0);
        
        Box blockDomainBox(Point::Ones(domainSize));
        err[nn] = 0;
        for (auto iter : layout)
        {
            int block = layout.block(iter);
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
                        // Apply Point
                        MBDataPoint dstDataPoint(iter, bi, layout);
                        MBPointInterpOp pointInterp(dstDataPoint, ghost[0], map, footprint, 4);
                        pointInterp.apply(hostDst, hostSrc);

                        double interpValue = hostDst[dstDataPoint](0);
                        double exactValue =  hostSrc[dstDataPoint](0);
                        double errorValue = interpValue - exactValue;
                        hostErr[dstDataPoint](0) = std::abs(errorValue);
                        err[nn] = std::max(std::abs(errorValue), err[nn]);
                    }
                }
            }
        }
#if PR_VERBOSE > 0
        std::cout << "error: " << err[nn] << std::endl;
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
        double rateErr = std::abs(rate - 4);
#if PR_VERBOSE > 0
        std::cout << "convergence rate: " << rate << std::endl;
#endif
        EXPECT_LT(rateErr, 0.1);
    }
}
#endif
#if 0
TEST(MBPointInterpOp, RigidXPoint) {
    pout() << "RIGID XPOINT TEST START" << std::endl;
    int domainSize = 8;
    int boxSize = 8;
    HDF5Handler h5;

    Array<double, DIM> exp{4,4,4,0,0,0};
    Array<double, DIM> offset{0,0,0,0,0,0};
    Array<double, DIM> k{1,1,1,1,1,1};

    // input footprint
    std::vector<Point> footprint;
    for (auto pi : Box::Kernel(3))
    {
        if (pi.abs().sum() <= 2)
        {
            footprint.push_back(pi);
        }
    }
    
    constexpr int N = 4;
    double err[N];
    for (int nn = 0; nn < N; nn++)
    {
        auto domain = buildXPoint(domainSize);
        Point boxSizeVect = Point::Ones(boxSize);
        MBDisjointBoxLayout layout(domain, boxSizeVect);

        Array<Point, DIM+1> ghost;
        ghost.fill(Point::Ones(4));
        ghost[0] = Point::Ones(1);
        Point boundGhost = Point::Ones(1);

        // initialize map
        MBMap<XPointMapRigid_t> map(XPointMapRigid, layout, ghost, boundGhost);

        // initialize data
        MBLevelBoxData<double, NCOMP, HOST> hostSrc(layout, ghost);
        MBLevelBoxData<double, NCOMP, HOST> hostDst(layout, ghost);
        MBLevelBoxData<double, NCOMP, HOST> hostErr(layout, ghost);

        hostSrc.initConvolve(f_polyM, map, exp, offset);
        hostSrc.fillBoundaries();
        hostDst.setVal(0);
        hostErr.setVal(0);
        
        Box blockDomainBox(Point::Ones(domainSize));
        err[nn] = 0;
        for (auto iter : layout)
        {
            int block = layout.block(iter);
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
                        // Apply Point
                        MBDataPoint dstDataPoint(iter, bi, layout);
                        MBPointInterpOp pointInterp(dstDataPoint, ghost[0], map, footprint, 4);
                        pointInterp.apply(hostDst, hostSrc);
                        
                        double interpValue = hostDst[dstDataPoint](0);
                        double exactValue =  hostSrc[dstDataPoint](0);
                        double errorValue = interpValue - exactValue;
                        hostErr[dstDataPoint](0) = std::abs(errorValue);
                        err[nn] = std::max(std::abs(errorValue), err[nn]);
                    }
                }
            }
        }
#if PR_VERBOSE > 0
        std::cout << "error: " << err[nn] << std::endl;
        h5.writeMBLevel({"err"}, map, hostErr, "MBInterpOpTests_XPointRigid_Err_%i", nn);
        h5.writeMBLevel({"phi"}, map, hostSrc, "MBInterpOpTests_XPointRigid_Src_%i", nn);
        h5.writeMBLevel({"phi"}, map, hostDst, "MBInterpOpTests_XPointRigid_Dst_%i", nn);
#endif
        domainSize *= 2;
        boxSize *= 2;
    }

    for (int ii = 1; ii < N; ii++)
    {
        double rate = log(err[ii-1]/err[ii])/log(2.0);
        double rateErr = std::abs(rate - 4);
#if PR_VERBOSE > 0
        std::cout << "convergence rate: " << rate << std::endl;
#endif
        EXPECT_LT(rateErr, 0.1);
    }
}
#endif
#if 0
TEST(MBPointInterpOp, XPointDisc) {
    pout() << "XPOINT DISC TEST START" << std::endl;
    int domainSize = 32;
    int boxSize = 32;
    HDF5Handler h5;

    Array<double, DIM> exp{4,4,4,0,0,0};
    Array<double, DIM> offset{0,0,0,0,0,0};
    Array<double, DIM> k{1,1,1,1,1,1};

    // input footprint
    std::vector<Point> footprint;
    for (auto pi : Box::Kernel(3))
    {
        if (pi.abs().sum() <= 2)
        {
            footprint.push_back(pi);
        }
    }
    
    constexpr int N = 3;
    double err[N];
    for (int nn = 0; nn < N; nn++)
    {
        auto domain = buildXPoint(domainSize);
        Point boxSizeVect = Point::Ones(boxSize);
        MBDisjointBoxLayout layout(domain, boxSizeVect);

        Array<Point, DIM+1> ghost;
        ghost.fill(Point::Ones(4));
        ghost[0] = Point::Ones(1);
        Point boundGhost = Point::Ones(1);

        // initialize map
        MBMap<XPointMap_t> map(XPointMap, layout, ghost, boundGhost);

        // initialize data
        MBLevelBoxData<double, NCOMP, HOST> hostSrc(layout, ghost);
        MBLevelBoxData<double, NCOMP, HOST> hostDst(layout, ghost);
        MBLevelBoxData<double, NCOMP, HOST> hostErr(layout, ghost);

        hostSrc.initConvolve(f_polyM, map, exp, offset);
        hostSrc.fillBoundaries();
        hostDst.setVal(0);
        hostErr.setVal(0);
        
        Box blockDomainBox(Point::Ones(domainSize));
        err[nn] = 0;
        for (auto iter : layout)
        {
            int block = layout.block(iter);
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
                        // Apply Point
                        MBDataPoint dstDataPoint(iter, bi, layout);
                        MBPointInterpOp pointInterp(dstDataPoint, ghost[0], map, footprint, 4);
                        pointInterp.apply(hostDst, hostSrc);
                        
                        double interpValue = hostDst[dstDataPoint](0);
                        double exactValue =  hostSrc[dstDataPoint](0);
                        double errorValue = interpValue - exactValue;
                        hostErr[dstDataPoint](0) = std::abs(errorValue);
                        err[nn] = std::max(std::abs(errorValue), err[nn]);
                    }
                }
            }
        }
#if PR_VERBOSE > 0
        std::cout << "error: " << err[nn] << std::endl;
        h5.writeMBLevel({"err"}, map, hostErr, "MBInterpOpTests_XPointDisc_Err_%i", nn);
        h5.writeMBLevel({"phi"}, map, hostSrc, "MBInterpOpTests_XPointDisc_Src_%i", nn);
        h5.writeMBLevel({"phi"}, map, hostDst, "MBInterpOpTests_XPointDisc_Dst_%i", nn);
#endif
        domainSize *= 2;
        boxSize *= 2;
    }

    for (int ii = 1; ii < N; ii++)
    {
        double rate = log(err[ii-1]/err[ii])/log(2.0);
        double rateErr = std::abs(rate - 4);
#if PR_VERBOSE > 0
        std::cout << "convergence rate: " << rate << std::endl;
#endif
        EXPECT_LT(rateErr, 0.1);
    }
}
#endif
#if 0
TEST(MBPointInterpOp, Ring) {
    pout() << "RING TEST START" << std::endl;
    int domainSize = 8;
    int boxSize = 8;
    HDF5Handler h5;

    Array<double, DIM> exp{4,4,4,0,0,0};
    Array<double, DIM> offset{0,0,0,0,0,0};
    Array<double, DIM> k{1,1,1,1,1,1};

    // input footprint
    std::vector<Point> footprint;
    for (auto pi : Box::Kernel(3))
    {
        if (pi.abs().sum() <= 2 && pi[1] == 0)
        //if (pi.abs().sum() <= 2)
        {
            footprint.push_back(pi);
        }
    }
#if PR_VERBOSE > 0
    pout() << "Input footprint: " << std::endl;
    for (auto pi : footprint)
    {
        pout() << pi << ", ";
    }
    pout() << std::endl;
#endif
    
    constexpr int N = 1;
    double errInf[N];
    double errL1[N];
    for (int nn = 0; nn < N; nn++)
    {
        Point domainSizeVect = Point::Ones(domainSize);
        Point boxSizeVect = Point::Ones(boxSize);
        for (int dir = 2; dir < DIM; dir++)
        {
            domainSizeVect[dir] = 2;
            boxSizeVect[dir] = 2;
        }
        auto domain = buildRing(domainSizeVect);
        MBDisjointBoxLayout layout(domain, boxSizeVect);

        Array<Point, DIM+1> ghost;
        ghost.fill(Point::Ones(4));
        ghost[0] = Point::Ones(1);
        Point boundGhost = Point::Ones(1);

        // initialize map
        MBMap<RingMap_t> map(RingMap, layout, ghost, boundGhost);

        // initialize data
        MBLevelBoxData<double, NCOMP, HOST> hostSrc(layout, ghost);
        MBLevelBoxData<double, NCOMP, HOST> hostDst(layout, ghost);
        MBLevelBoxData<double, NCOMP, HOST> hostErr(layout, ghost);

        hostSrc.initConvolve(f_polyM, map, exp, offset);
        hostSrc.fillBoundaries();
        hostDst.setVal(0);
        hostErr.setVal(0);
        
        Box blockDomainBox(domainSizeVect);
        errInf[nn] = 0;
        errL1[nn] = 0;
        for (auto iter : layout)
        {
            int block = layout.block(iter);
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

                        // Apply Point
                        MBDataPoint dstDataPoint(iter, bi, layout);
                        MBPointInterpOp pointInterp(dstDataPoint, ghost[0], map, footprint, 4);
                        pointInterp.apply(hostDst, hostSrc);

                        pout() << "Coefs at point " << bi << std::endl;
                        auto coefs = pointInterp.coefs(hostSrc);
                        coefs.print("%10.2e");

                        double interpValue = hostDst[dstDataPoint](0);
                        double exactValue =  hostSrc[dstDataPoint](0);
                        double errorValue = interpValue - exactValue;
                        hostErr[dstDataPoint](0) = std::abs(errorValue);
                        errInf[nn] = std::max(std::abs(errorValue), errInf[nn]);
                        errL1[nn] += std::abs(errorValue);
                    }
                }
            }
        }
        double surfaceArea = std::pow(domainSize, DIM-1);
        errL1[nn] /= surfaceArea;
#if PR_VERBOSE > 0
        std::cout << "error (inf): " << errInf[nn] << std::endl;
        std::cout << "error  (L1): " << errL1[nn] << std::endl;
        h5.writeMBLevel({"err"}, map, hostErr, "MBInterpOpTests_Ring_Err_%i", nn);
        h5.writeMBLevel({"phi"}, map, hostSrc, "MBInterpOpTests_Ring_Src_%i", nn);
        h5.writeMBLevel({"phi"}, map, hostDst, "MBInterpOpTests_Ring_Dst_%i", nn);
        h5.writeMBLevel({"J"}, map, map.jacobian(), "MBInterpOpTests_Ring_J_%i", nn);
#endif
        domainSize *= 2;
        boxSize *= 2;
    }

    for (int ii = 1; ii < N; ii++)
    {
        double rateInf = log(errInf[ii-1]/errInf[ii])/log(2.0);
        double rateInfErr = std::abs(rateInf - 4);
        double rateL1 = log(errL1[ii-1]/errL1[ii])/log(2.0);
        double rateL1Err = std::abs(rateL1 - 4);
#if PR_VERBOSE > 0
        std::cout << "convergence rate (inf): " << rateInf << std::endl;
        std::cout << "convergence rate (L1): " << rateL1 << std::endl;
#endif
        EXPECT_LT(rateInfErr, 0.1);
        EXPECT_LT(rateL1Err, 0.1);
    }
}
#endif
#if DIM == 3
#if 1
TEST(MBPointInterpOp, SphericalShell) {
    int domainSize = 8;
    int boxSize = 8;
    int thickness = 1;
    HDF5Handler h5;

    Array<double, DIM> exp{0,0,0,0,0,0};
    Array<double, DIM> offset{0,0,0,0,0,0};
    Array<double, DIM> k{1,1,1,1,1,1};

    // input footprint
    std::vector<Point> footprint;
    for (auto pi : Box::Kernel(2))
    {
        if (pi.abs().sum() <= 2)// && (pi[0] == 0))
        {
            footprint.push_back(pi);
        }
    }
    std::vector<Point> exponents;
    std::vector<std::string> coefVarNames;
    for (auto bi : Box::Cube(4))
    {
        if (bi.sum() < 4)
        {
            exponents.push_back(bi);
            std::string varname;
            std::string xyz = "xyz";
            for (int dir = 0; dir < DIM; dir++)
            {
                if (bi[dir] == 0) {continue;}
                else {
                    if (varname.size() > 0) {varname+="*";}
                    varname += xyz[dir];
                    if (bi[dir] > 1)
                    {
                        varname+=("^" + std::to_string(bi[dir]));
                    }
                }
            }
            if (varname.size() == 0) {varname = "1";}
            coefVarNames.push_back(varname);
        }
    }
   
#if PR_VERBOSE > 0
    pout() << "Input footprint: " << std::endl;
    for (auto pi : footprint)
    {
        pout() << pi << ", ";
    }
    pout() << std::endl;
#endif
    constexpr int N = 3;
    double err[N];
    for (int nn = 0; nn < N; nn++)
    {
        //auto domain = buildPartialThinCubeSphere(domainSize, thickness);
        auto domain = buildThinCubeSphere(domainSize, thickness);
        Point boxSizeVect = Point::Ones(boxSize);
        boxSizeVect[2] = thickness;
        MBDisjointBoxLayout layout(domain, boxSizeVect);

        Array<Point, DIM+1> ghost;
        ghost.fill(Point::Ones(4));
        ghost[0] = Point::Ones(1);
        Point boundGhost = Point::Ones(1);
        ghost[0][2] = 0;

        // initialize map
        //MBMap<PartialThinCubedSphereMap_t> map(PartialThinCubedSphereMap, layout, ghost, boundGhost);
        MBMap<ThinCubedSphereMap_t> map(ThinCubedSphereMap, layout, ghost, boundGhost);
        
        // initialize data
        MBLevelBoxData<double, NCOMP, HOST> hostSrc(layout, ghost);
        MBLevelBoxData<double, NCOMP, HOST> hostDst(layout, ghost);
        MBLevelBoxData<double, NCOMP, HOST> hostErr(layout, ghost);
        MBLevelBoxData<double, 20, HOST> hostCoefs(layout, ghost);

        hostSrc.initConvolve(f_polyM, map, exp, offset);
        hostSrc.fillBoundaries();
        hostDst.setVal(0);
        hostErr.setVal(0);
        hostCoefs.setVal(0);
        
        Box blockDomainBox(Point::Ones(domainSize));
        err[nn] = 0;
        for (auto iter : layout)
        {
            int block = layout.block(iter);
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
                        // Apply Point
                        MBDataPoint dstDataPoint(iter, bi, layout);
                        MBPointInterpOp pointInterp(dstDataPoint, ghost[0], map, footprint, 4);
                        pointInterp.apply(hostDst, hostSrc);

                        auto coefs = pointInterp.coefs(hostSrc);
                        for (int ei = 0; ei < 20; ei++)
                        {
                            hostCoefs[iter](bi,ei) = coefs(ei,0);
                        }
                        
                        double interpValue = hostDst[dstDataPoint](0);
                        double exactValue =  hostSrc[dstDataPoint](0);
                        double errorValue = interpValue - exactValue;
                        hostErr[dstDataPoint](0) = std::abs(errorValue);
                        err[nn] = std::max(std::abs(errorValue), err[nn]);
                    }
                }
            }
        }
#if PR_VERBOSE > 0
        std::cout << "error: " << err[nn] << std::endl;
        h5.writeMBLevel({"err"}, map, hostErr, "MBInterpOpTests_SphericalShell_Err_%i", nn);
        h5.writeMBLevel({"phi"}, map, hostSrc, "MBInterpOpTests_SphericalShell_Src_%i", nn);
        h5.writeMBLevel({"phi"}, map, hostDst, "MBInterpOpTests_SphericalShell_Dst_%i", nn);
        h5.writeMBLevel({"J"}, map, map.jacobian(), "MBInterpOpTests_SphericalShell_J_%i", nn);
        h5.writeMBLevel(coefVarNames, map, hostCoefs, "MBInterpOpTests_SphericalShell_Coefs_%i", nn);
#endif
        domainSize *= 2;
        boxSize *= 2;
    }
    for (int ii = 1; ii < N-1; ii++)
    {
        for (int mi = 0; mi < 20; mi++)
        {
#if PR_VERBOSE > 0
            //std::cout << "Error in moment " << exponents[mi] << ": " << momentErrNorm[ii-1][mi] << std::endl;
#endif
        }
    }

    for (int ii = 1; ii < N; ii++)
    {
        double rate = log(err[ii-1]/err[ii])/log(2.0);
        double rateErr = std::abs(rate - 4);
#if PR_VERBOSE > 0
        std::cout << "convergence rate: " << rate << std::endl;
#endif
        EXPECT_LT(rateErr, 0.1);
    }
}
#endif
#if 0
TEST(MBPointInterpOp, PolarShell) {
    pout() << "PolarShell" << std::endl;
    int domainSize = 8;
    int boxSize = 8;
    int thickness = 1;
    HDF5Handler h5;

    Array<double, DIM> exp{4,4,4,0,0,0};
    Array<double, DIM> offset{0.1,0.1,0.1,0,0,0};
    Array<double, DIM> k{1,1,1,1,1,1};

    // input footprint
    std::vector<Point> footprint;
    for (auto pi : Box::Kernel(2))
    {
        if (pi.abs().sum() <= 2 && (pi[0] == 0))
        {
            footprint.push_back(pi);
        }
    }
   
#if PR_VERBOSE > 0
    pout() << "Input footprint: " << std::endl;
    for (auto pi : footprint)
    {
        pout() << pi << ", ";
    }
    pout() << std::endl;
#endif
    constexpr int N = 2;
    double err[N];
    for (int nn = 0; nn < N; nn++)
    {
        auto domain = buildPolarShell(domainSize, thickness);
        Point boxSizeVect = Point::Ones(boxSize);
        boxSizeVect[0] = thickness;
        MBDisjointBoxLayout layout(domain, boxSizeVect);

        Array<Point, DIM+1> ghost;
        ghost.fill(Point::Ones(4));
        ghost[0] = Point::Ones(1);
        Point boundGhost = Point::Ones(1);

        // initialize map
        MBMap<MHDMap_t> map(MHDMap, layout, ghost, boundGhost);

        // initialize data
        MBLevelBoxData<double, NCOMP, HOST> hostSrc(layout, ghost);
        MBLevelBoxData<double, NCOMP, HOST> hostDst(layout, ghost);
        MBLevelBoxData<double, NCOMP, HOST> hostErr(layout, ghost);

        //hostSrc.initConvolve(f_polyM, map, exp, offset);
        hostSrc.initialize(f_phiM, map, k, offset);
        hostSrc.fillBoundaries();
        hostDst.setVal(0);
        hostErr.setVal(0);
        
        Box blockDomainBox(Point::Ones(domainSize));
        err[nn] = 0;
        for (auto iter : layout)
        {
            int block = layout.block(iter);
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
                        // Apply Point
                        MBDataPoint dstDataPoint(iter, bi, layout);
                        MBPointInterpOp pointInterp(dstDataPoint, ghost[0], map, footprint, 4);
                        pointInterp.apply(hostDst, hostSrc);

                        //pout() << "Coefs at point " << bi << std::endl;
                        //auto coefs = pointInterp.coefs(hostSrc);
                        //coefs.print("%10.2e");

                        double interpValue = hostDst[dstDataPoint](0);
                        double exactValue =  hostSrc[dstDataPoint](0);
                        double errorValue = interpValue - exactValue;
                        hostErr[dstDataPoint](0) = std::abs(errorValue);
                        err[nn] = std::max(std::abs(errorValue), err[nn]);
                    }
                }
            }
        }
#if PR_VERBOSE > 0
        std::cout << "error: " << err[nn] << std::endl;
        h5.writeMBLevel({"err"}, map, hostErr, "MBInterpOpTests_PolarShell_Err_%i", nn);
        h5.writeMBLevel({"phi"}, map, hostSrc, "MBInterpOpTests_PolarShell_Src_%i", nn);
        h5.writeMBLevel({"phi"}, map, hostDst, "MBInterpOpTests_PolarShell_Dst_%i", nn);
        h5.writeMBLevel({"J"}, map, map.jacobian(), "MBInterpOpTests_PolarShell_J_%i", nn);
#endif
        domainSize *= 2;
        boxSize *= 2;
        thickness *= 2;
    }

    for (int ii = 1; ii < N; ii++)
    {
        double rate = log(err[ii-1]/err[ii])/log(2.0);
        double rateErr = std::abs(rate - 4);
#if PR_VERBOSE > 0
        std::cout << "convergence rate: " << rate << std::endl;
#endif
        EXPECT_LT(rateErr, 0.1);
    }
}
#endif
#endif

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
        }
    }
 
    M[1].define(1, exponents.size());
    for (auto ei = 0; ei < exponents.size(); ei++)
    {
        M[1](0, ei) = Xp_avg[0][ei](a_dst);
    }

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
    for (auto pi : Box::Kernel(2))
    {
        if (pi.abs().sum() <= 2)
        {
            footprint.push_back(pi);
            std::cout << pi << std::endl;
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
        MBDataPoint t(*layout.begin(), pi, layout);
        auto& op = interp(t);
        auto M = computeM(map, boxSize, pi);
        auto Cinv = M[0].inverse();
        auto S = M[1]*Cinv;
        auto EC = M[0] - op.MC();
        auto ED = M[1] - op.MD();
        auto ES = S - op.MS();
        //EC.print();
        //ED.print();
        //ES.print();
        CErr = std::max(CErr, (EC.absMax()));
        DErr = std::max(DErr, (ED.absMax()));
        SErr = std::max(SErr, (ES.absMax()));
        //M[0].print();
        //M[1].print();
        pout() << "-----------------------------------------------------" << std::endl;
        pout() << "Stencil at point: " << pi << std::endl;
        pout() << "Source data: " << std::endl;
        pout() << "Sum of stencil coefficients: " << op.MS().sum() << std::endl;
        for (auto& si : op.sources())
        {
            pout() << si.point << ", ";
        }
        pout() << std::endl;
        S.print();
    }
    
    EXPECT_LT(CErr, 1e-12);
    EXPECT_LT(DErr, 1e-12);
    EXPECT_LT(SErr, 1e-12);

#if PR_VERBOSE > 0
    std::cout << "Error in C: " << CErr << std::endl;;
    std::cout << "Error in D: " << DErr << std::endl;;
    std::cout << "Error in S: " << SErr << std::endl;;
#endif
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
