#include <gtest/gtest.h>
#include "ProtoAMR.H"
#include "ProtoMMB.H"
#include "Lambdas.H"

using namespace Proto;

DisjointBoxLayout testLayout(int domainSize, Point boxSize)
{
    Box domainBox = Box::Cube(domainSize); 
    Box patchBox = domainBox.coarsen(boxSize);
    std::vector<Point> patches;
    for (auto patch : patchBox)
    {
        patches.push_back(patch);
        //if (patch != Point::Zeros()) { patches.push_back(patch); }
    }
    std::array<bool, DIM> periodicity;
    periodicity.fill(true);
    ProblemDomain domain(domainBox, periodicity);
    return DisjointBoxLayout(domain, patches, boxSize);
}

template<MemType MEM>
PROTO_KERNEL_START
void f_testMapF(Point& a_pt, Var<double,3,MEM>& a_X)
{
    
    double x = a_pt[0];
    double y = a_pt[1] + x;
    a_X(0) = x;
    a_X(1) = y;
    a_X(2) = 0;
#if DIM > 2
    a_X(2) = a_pt[2];
#endif

}
PROTO_KERNEL_END(f_testMapF, f_testMap);

TEST(HDF5, MMB) {
    int domainSize = 64;
    Point boxSize = Point::Ones(domainSize);
    double dx = 1.0/domainSize;
    auto layout = testLayout(domainSize, boxSize);
    LevelBoxData<double, 1, HOST> data(layout, Point::Zeros());
    LevelBoxData<double, 3, HOST, PR_NODE> map(layout, Point::Zeros());
    data.initialize(f_phi, dx);
    map.initialize(f_testMap);
    HDF5Handler h5;
    h5.writeLevel(dx, data, "DATA");
    h5.writeLevel({"x", "y", "z"}, dx, map,"DATA.map");
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
