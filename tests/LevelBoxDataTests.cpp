#include <gtest/gtest.h>
#include "Proto.H"

using namespace Proto;

TEST(LevelBoxData, Initialize) {
    int domainSize = 64;
    int boxSize = 16;
    double dx = 1.0/domainSize;
    double offset = 0.125;
    Point boxSizeVect = Point::Ones(boxSize);
    Box domainBox = Box::Cube(domainSize);
    Box patchBox = domainBox.coarsen(boxSizeVect);
    std::vector<Point> patches;
    for (auto patch : patchBox)
    {
        if (patch[0] != 0) { patches.push_back(patch); }
    }
    std::array<bool, DIM> periodicity;
    periodicity.fill(true);
    ProblemDomain domain(domainBox, periodicity);
    DisjointBoxLayout layout(domain, patches, boxSizeVect);
    LevelBoxData<double, 1, HOST> hostData(layout, Point::Ones());
    hostData.initialize(f_phi, dx, offset);
#ifdef PROTO_CUDA
    LevelBoxData<double, 1, DEVICE> deviData(layout, Point::Ones());
    deviData.initialize(f_phi, dx, offset);
#endif
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
