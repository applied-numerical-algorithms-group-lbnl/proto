#include <gtest/gtest.h>
#include "Proto.H"
#include "Lambdas.H"

#define NCOMP 1

using namespace Proto;

TEST(MBPointInterpOp, XPointApply) {
    int domainSize = 64;
    int boxSize = 32;
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

        Point k{1,2,3,4,5,6};
        Array<double, DIM> offset{1,1,1,1,1,1};
        hostSrc.initialize(f_phiM, map, k, offset);
        hostSrc.fillBoundaries();
        hostDst.setVal(0);
        hostErr.setVal(0);

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
                if (dir != Point::Basis(0) && dir != Point::Basis(1)
                        && dir != (Point::Basis(1) + Point::Basis(0))){continue;}
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
        h5.writeMBLevel({"err"}, map, hostErr, "Err_%i", nn);
        pout() << "Error: " << err[nn] << std::endl;
        domainSize *= 2;
        boxSize *= 2;
    }

    double rateTol = 0.1;
    for (int ii = 1; ii < N; ii++)
    {
        double rate = log(err[ii-1]/err[ii])/log(2.0);
        double rateErr = abs(rate - 4);
        EXPECT_LT(rateErr, rateTol);
        pout() << "Convergence Rate: " << log(err[ii-1]/err[ii])/log(2.0) << std::endl;
    }
}

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
