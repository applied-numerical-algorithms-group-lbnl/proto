#include <gtest/gtest.h>
#include "Proto.H"
#include "Lambdas.H"

#define NCOMP 1

using namespace Proto;



TEST(MBPointInterpOp, Constructor) {
    HDF5Handler h5;

    int domainSize = 32;
    int boxSize = 16;
#if DIM > 2
    domainSize = 8;
    boxSize = 4;
#endif
    constexpr int N = 6;
    Array<double, N> edgeErr(0);
    Array<double, N> cornErr(0);
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
        MBMap map(XPointMapRigid, layout, ghost, boundGhost);

        // initialize data
        MBLevelBoxData<double, NCOMP, HOST> hostSrc(layout, ghost, boundGhost);
        MBLevelBoxData<double, NCOMP, HOST> hostDst(layout, ghost, boundGhost);
        MBLevelBoxData<double, NCOMP, HOST> hostErr(layout, ghost, boundGhost);


        h5.writeMBLevel({"x", "y", "z"}, map.map(), "J_N%i.map",nn);
        h5.writeMBLevel({"J"}, map.jacobian(), "J_N%i",nn);
        

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

        // inputs
        Point p0 = Point::Ones(domainSize-1);// + Point::Basis(0);
        Point patchID = Point::Ones((domainSize / boxSize) - 1);
        auto mbIndex = layout.find(patchID, 0);

        MBDataPoint dstPoint(mbIndex, p0);

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
                auto patchPoint = layout.point(iter);
                if (patchPoint == (Point::Basis(0) + Point::Basis(1)))
                {
                    cornErr[nn] = e > cornErr[nn] ? e : cornErr[nn];
                } else {
                    edgeErr[nn] = e > edgeErr[nn] ? e : edgeErr[nn];
                }
            }
        }
        std::cout << "Edge Error: " << edgeErr[nn] << std::endl;
        std::cout << "Corner Error: " << cornErr[nn] << std::endl;
        h5.writeMBLevel({"x", "y", "z"}, map.map(), "INTERP_N%i.map",nn);
        h5.writeMBLevel({"err"}, hostErr, "INTERP_N%i",nn);
        domainSize *= 2;
        boxSize *= 2;
    }

    for (int ii = 1; ii < N; ii++)
    {
        std::cout << "Edge Rate: " << log(edgeErr[ii-1]/edgeErr[ii])/log(2.0) << std::endl;
        std::cout << "Corner Rate: " << log(cornErr[ii-1]/cornErr[ii])/log(2.0) << std::endl;
    }
}

int main(int argc, char *argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
#ifdef PR_MPI
    MPI_Init(&argc, &argv);
#endif
    PR_TIMER_SETFILE("DIM=" + to_string(DIM) + ".numProc=" + to_string(numProc())
        + "MBInterpOpTests.time.table");
    int result = RUN_ALL_TESTS();
    PR_TIMER_REPORT();
#ifdef PR_MPI
    MPI_Finalize();
#endif
    return result;
}
