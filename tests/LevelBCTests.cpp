#include <gtest/gtest.h>
#include "Proto.H"
#include "Lambdas.H"
#include "LevelBC_Constant.H"

using namespace Proto;

TEST(LevelBC, Constant) {
    HDF5Handler h5;
    int domainSize = 16;
    int boxSize = 8;

    Box domainBox = Box::Cube(domainSize); 
    std::array<bool, DIM> periodicity;
    periodicity.fill(false);
    ProblemDomain domain(domainBox, periodicity);
    DisjointBoxLayout layout(domain, Point::Ones(boxSize));

    LevelBoxData<double, 1, HOST> data(layout, Point::Ones());
    data.setVal(0);
    
    ConstantBC<double, 1, HOST, PR_CELL> bc(layout);
    bc.setNumConstants(1);
    bc.setConstant(0, 7);
#if PR_VERBOSE > 0
    h5.writeLevel(1, data, "LevelBCTests_Data_0");
#endif
    bc.apply(data);
#if PR_VERBOSE > 0
    h5.writeLevel(1, data, "LevelBCTests_Data_1");
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
