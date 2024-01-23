#include <gtest/gtest.h>
#include "Proto.H"
#include "Lambdas.H"
#include "LevelBC_Constant.H"

using namespace Proto;

#if 0
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
#endif
TEST(LevelBC, ConstDirichlet)
{
    HDF5Handler h5;
    int domainSize = 16;
    int boxSize = 8;
    int ghostSize = 1;
    double bcValue = 3.0;

    Box domainBox = Box::Cube(domainSize); 
    std::array<bool, DIM> periodicity;
    periodicity.fill(false);
#if DIM > 1
    periodicity[1] = true;
#endif
    ProblemDomain domain(domainBox, periodicity);
    DisjointBoxLayout layout(domain, Point::Ones(boxSize));

    LevelBoxData<double, DIM, HOST> data_0(layout, Point::Ones(ghostSize));
    LevelBoxData<double, DIM, HOST> data_1(layout, Point::Ones(ghostSize));
    data_0.setVal(7);
    data_1.setVal(7);

    std::cout << "DIM: " << DIM << std::endl;
    ConstDirichletBC<double, DIM, HOST, PR_CELL> globalBC(layout);
    globalBC.setVal(bcValue);
    ConstDirichletBC<double, DIM, HOST, PR_CELL> localBC(layout);
    localBC.setVal(bcValue, 0);

#if PR_VERBOSE > 0
    h5.writeLevel(1, data_0, "LevelBCTests_GlobalData_0");
    h5.writeLevel(1, data_1, "LevelBCTests_LocalData_0");
#endif
    globalBC.apply(data_0);
    localBC.apply(data_1);
#if PR_VERBOSE > 0
    h5.writeLevel(1, data_0, "LevelBCTests_GlobalData_1");
    h5.writeLevel(1, data_1, "LevelBCTests_LocalData_1");
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
