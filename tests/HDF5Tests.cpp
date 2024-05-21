#include <gtest/gtest.h>
#include "Proto.H"
#include "Lambdas.H"

#define XPOINT_SIZE 5
#define NCOMP 2

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

TEST(HDF5, ReadMBLevel)
{
    HDF5Handler h5;
    int domainSize = 64;
    int boxSize = 16;
    auto domain = buildXPoint(domainSize);
    Point boxSizeVect = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);
    MBLevelBoxData<double, 3, HOST> data(layout, Point::Ones());
    data.setVal(7);
    h5.writeMBLevel({"var1", "var2","var3"}, data, "DATA_0");
    
    ProblemDomain readDomain(Box(Point(64,64)), false);
    DisjointBoxLayout readLayout(readDomain, Point::Ones(16));

    LevelBoxData<double, 3, HOST> readData(readLayout, Point::Ones());
    readData.setVal(-7);
    h5.writeLevel({"var1", "var2","var3"}, readData, "DATA_1");
    h5.readLevel(readData, "DATA_0");

    h5.writeLevel({"var1", "var2","var3"}, readData, "DATA_2");
}

/*
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
    map.initialize(f_XPointMap, XPOINT_SIZE, domainSize);

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

    HDF5Handler h5;
    int domainSize = 64;
    int boxSize = 16;
    auto domain = buildCubeSphere(domainSize);
    Point boxSizeVect = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);
    std::array<Point, DIM+1> ghost;
    ghost.fill(Point::Ones());
    MBLevelBoxData<double, NCOMP, HOST> data(layout, Point::Ones());
    data.initialize(f_MBPointID);
    MBLevelBoxData<double, 3, HOST, PR_NODE> map(layout, Point::Ones());
    map.initialize(f_CubeSphereMap, Point::Ones(domainSize), 1.0, 2.0);
    h5.writeMBLevel({"var1", "var2"}, data, "CUBE_SPHERE");
    h5.writeMBLevel({"x", "y", "z"}, map, "CUBE_SPHERE.map");
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
