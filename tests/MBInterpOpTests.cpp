#include <gtest/gtest.h>
#include "Proto.H"
#include "Lambdas.H"

#define NCOMP 1

using namespace Proto;



TEST(MBPointInterpOp, Constructor) {
    HDF5Handler h5;

    int domainSize = 16;
    int boxSize = 8;
    auto domain = buildXPoint(domainSize);
    Point boxSizeVect = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);

    Array<Point, DIM+1> ghost;
    ghost.fill(Point::Ones(5));
    ghost[0] = Point::Ones(2);
    Point boundGhost = Point::Ones();
    // initialize data
    MBLevelBoxData<double, NCOMP, HOST> hostSrc(layout, ghost, boundGhost);
    MBLevelBoxData<double, NCOMP, HOST> hostDst(layout, ghost, boundGhost);
    hostSrc.initialize(f_MBPointID);
    hostSrc.fillBoundaries();
    hostDst.setVal(0);
    
    // initialize map
    MBMap map(XPointMap, layout, ghost, boundGhost);

        

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

#if 0
    MBPointInterpOp interp(dstPoint, map, footprint, 4);
    interp.apply(hostDst, hostSrc); 
#else
    MBInterpOp interp(map, footprint, ghost[0], 4);
    interp.apply(hostDst, hostSrc);
#endif
    h5.writeMBLevel({"x", "y", "z"}, map.map(), "INTERP.map");
    h5.writeMBLevel({"x", "y", "z"}, hostDst, "INTERP");

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
