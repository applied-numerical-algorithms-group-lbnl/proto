#include <gtest/gtest.h>
#include "Proto.H"
#include "Lambdas.H"

using namespace Proto;

TEST(FluxRegister, Construction) {
    constexpr unsigned int C = 2;
    int domainSize = 32;
    int numLevels = 2;
    double offset = 0.125;
    double dx = 1.0/domainSize;
    Point refRatio = Point::Ones(2);
    Point boxSize = Point::Ones(16);
    auto grid = telescopingGrid(domainSize, numLevels, refRatio, boxSize);
    
    std::array<double, DIM> dxVect;
    dxVect.fill(dx);
    LevelFluxRegister<double, C, HOST> fr(grid[0], grid[1], refRatio, dxVect);
}

TEST(FluxRegister, RefluxCenter) {
    constexpr unsigned int C = 2;
    int domainSize = 32;
    int numLevels = 2;
    int ghostSize = 2;
    double offset = 0.125;
    double dx = 1.0/domainSize;
    Point refRatio = Point::Ones(2);
    Point boxSize = Point::Ones(16);
    auto grid = telescopingGrid(domainSize, numLevels, refRatio, boxSize);

    LevelBoxData<double, C, HOST> L0(grid[0], Point::Ones(ghostSize));
    LevelBoxData<double, C, HOST> L1(grid[1], Point::Ones(ghostSize));

    L0.setVal(0);
    L1.setVal(0);

    std::array<double, DIM> dxVect;
    dxVect.fill(dx);
    LevelFluxRegister<double, C, HOST> fr(grid[0], grid[1], refRatio, dxVect);

    for (auto citer : grid[0])
    {
        for (int dir = 0; dir < DIM; dir++)
        {
            BoxData<double, C, HOST> flux(grid[0][citer].extrude(Point::Basis(dir,1)));
            flux.setVal(dir+1);
            fr.incrementCoarse(flux, citer, 1.0, dir);
        }
    }
    
    for (auto fiter : grid[1])
    {
        for (int dir = 0; dir < DIM; dir++)
        {
            BoxData<double, C, HOST> flux(grid[1][fiter].extrude(Point::Basis(dir,1)));
            flux.setVal((dir+1)*10);
            fr.incrementFine(flux, fiter, 1.0, dir);
        }
    }

#if PR_VERBOSE > 0
    HDF5Handler h5;
    h5.writeLevel(dx, L0, "LevelFluxRegisterTests_L0_0");
#endif
    fr.reflux(L0, 1.0);
#if PR_VERBOSE > 0
    h5.writeLevel(dx, L0, "LevelFluxRegisterTests_L0_1");
#endif
        

}

TEST(FluxRegister, RefluxCorner) {
    constexpr unsigned int C = 2;
    int domainSize = 32;
    int numLevels = 2;
    int ghostSize = 2;
    double offset = 0.125;
    double dx = 1.0/domainSize;
    Point refRatio = Point::Ones(2);
    Point boxSize = Point::Ones(16);
    auto grid = cornerGrid(domainSize, numLevels, refRatio, boxSize);

    LevelBoxData<double, C, HOST> L0(grid[0], Point::Ones(ghostSize));
    LevelBoxData<double, C, HOST> L1(grid[1], Point::Ones(ghostSize));

    L0.setVal(0);
    L1.setVal(0);

    std::array<double, DIM> dxVect;
    dxVect.fill(dx);
    LevelFluxRegister<double, C, HOST> fr(grid[0], grid[1], refRatio, dxVect);

    for (auto citer : grid[0])
    {
        for (int dir = 0; dir < DIM; dir++)
        {
            BoxData<double, C, HOST> flux(grid[0][citer].extrude(Point::Basis(dir,1)));
            flux.setVal(dir+1);
            fr.incrementCoarse(flux, citer, 1.0, dir);
        }
    }
    
    for (auto fiter : grid[1])
    {
        for (int dir = 0; dir < DIM; dir++)
        {
            BoxData<double, C, HOST> flux(grid[1][fiter].extrude(Point::Basis(dir,1)));
            flux.setVal((dir+1)*10);
            fr.incrementFine(flux, fiter, 1.0, dir);
        }
    }

#if PR_VERBOSE > 0
    HDF5Handler h5;
    h5.writeLevel(dx, L0, "LevelFluxRegisterTests_L0_0");
    h5.writeLevel(dx/refRatio[0], L1, "LevelFluxRegisterTests_L1_0");
#endif
    fr.reflux(L0, 1.0);
#if PR_VERBOSE > 0
    h5.writeLevel(dx, L0, "LevelFluxRegisterTests_L0_1");
    h5.writeLevel(dx/refRatio[0], L1, "LevelFluxRegisterTests_L1_1");
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
