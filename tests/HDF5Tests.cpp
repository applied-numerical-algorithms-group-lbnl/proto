#include <gtest/gtest.h>
#define XPOINT_SIZE 5
#define NCOMP 2
#include "ProtoMMB.H"
#include "ProtoAMR.H"
#include "Lambdas.H"


using namespace Proto;

MBProblemDomain buildXPoint(int a_domainSize)
{
    MBProblemDomain domain(XPOINT_SIZE);
    auto CCW = CoordPermutation::ccw();
    for (int ii = 0; ii < XPOINT_SIZE; ii++)
    {
        domain.defineBoundary(ii, (ii+1) % XPOINT_SIZE, 0, Side::Hi, CCW);
    }
    for (int bi = 0; bi < XPOINT_SIZE; bi++)
    {
        domain.defineDomain(bi, Point::Ones(a_domainSize));
    }
    return domain;
}

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

TEST(HDF5, MMBOffsets)
{
    HDF5Handler h5;
    int domainSize = 64;
    int boxSize = 16;
    auto domain = buildXPoint(domainSize);
    Point boxSizeVect = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);
    std::array<Point, DIM+1> ghost;
    ghost.fill(Point::Ones());
    MBLevelBoxData<double, NCOMP, HOST> data(layout, Point::Ones());
    data.initialize(f_MBPointID);
    MBLevelBoxData<double, 3, HOST, PR_NODE> map(layout, Point::Ones());
    map.initialize(f_XPointMap, domainSize);

    for (int bi = 0; bi < XPOINT_SIZE; bi++)
    {
        data.blockData(bi).setVal(bi,0);
        if (NCOMP > 1)
        {
            data.blockData(bi).setVal(bi+10,1);
        }
   //     h5.writeLevel(1.0, data.blockData(bi), "DATA_%i", bi);
    }
    
    unsigned int numBoxes = pow(domainSize / boxSize, DIM)*XPOINT_SIZE;
    std::vector<unsigned int> boxesPerProc(numProc(), numBoxes / numProc());
    for (int ii = 0; ii < numBoxes % numProc(); ii++) { boxesPerProc[ii]++; }
    
    size_t layoutOffset = 0;
    size_t dataOffset = 0;
    for (int ii = 0; ii < numProc(); ii++)
    {
        EXPECT_EQ(layoutOffset, layout.offset(ii));
        EXPECT_EQ(dataOffset, data.offset(ii));
        layoutOffset += boxesPerProc[ii];
        dataOffset += boxesPerProc[ii]*data.blockData(0).patchSize();
    }
    h5.writeMBLevel({"var1", "var2"}, data, "DATA");
    h5.writeMBLevel({"x", "y", "z"}, map, "DATA.map");
}
TEST(HDF5, CubeSphere)
{
    int domainSize = 64;
    Point boxSize = Point::Ones(domainSize);
    double dx = 1.0/domainSize;
    auto layout = testLayout(domainSize, boxSize);
    LevelBoxData<double, 1, HOST> data(layout, Point::Zeros());
    LevelBoxData<double, 3, HOST, PR_NODE> map(layout, Point::Zeros());
    data.initialize(f_phi, dx);
    map.initialize(f_CubeSphereMap, 0, boxSize, 1.0);
    HDF5Handler h5;
    h5.writeLevel(dx, data, "CUBE_SPHERE");
    h5.writeLevel({"x", "y", "z"}, dx, map,"CUBE_SPHERE.map");

}
/*
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
*/

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
