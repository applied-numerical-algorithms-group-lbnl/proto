#include <gtest/gtest.h>
#include "Proto.H"
#include "Lambdas.H"
#include "MBLevelMap_Shear.H"
#include "MBLevelMap_CubeSphereShell.H"

using namespace Proto;

TEST(MBInterpOp, ShearTest)
{
    int domainSize = 8;
    int boxSize = 8;
    Array<double, DIM> exp{0,0,0,0,0,0};
    Array<double, DIM> offset{0,0,0,0,0,0};
    HDF5Handler h5;

    auto domain = buildShear(domainSize);
    Point boxSizeVect = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);

    Array<Point, DIM+1> ghost;
    ghost.fill(Point::Ones(4));
    ghost[0] = Point::Ones(1);

    // initialize map
    MBLevelMap_Shear<HOST> map;
    map.define(layout, ghost);

    MBLevelBoxData<double, 1, HOST> hostSrc(layout, ghost);
    hostSrc.initConvolve(f_polyM, map, exp, offset);

    // input footprint
    std::vector<Point> footprint;
    for (auto pi : Box::Kernel(3))
    {
        if (pi.abs().sum() <= 2)
        {
            footprint.push_back(pi);
        }
    }
    
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
                }
            }
        }
    }
}

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
