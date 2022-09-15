#include <gtest/gtest.h>
#include "Proto.H"
#include "Lambdas.H"

using namespace Proto;

TEST(FluxRegister, Construction) {
    constexpr unsigned int C = 2;
    HDF5Handler h5;
    int domainSize = 32;
    int numLevels = 2;
    double offset = 0.125;
    double dx = 1.0/domainSize;
    Point refRatio = Point::Ones(2);
    Point boxSize = Point::Ones(16);
    auto grid = telescopingGrid(domainSize, numLevels, refRatio, boxSize);
    
    std::array<double, DIM> dxVect;
    dxVect.fill(dx);
    //LevelFluxRegister<double, C, HOST>
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
