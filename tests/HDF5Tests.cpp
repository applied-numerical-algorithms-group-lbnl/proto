#include <gtest/gtest.h>
#include "Proto.H"
#include "TestFunctions.H"

#define XPOINT_SIZE 5
#define NCOMP 2

using namespace Proto;

namespace {
    DisjointBoxLayout testLayout(int domainSize, Point boxSize)
    {
        Box domainBox = Box::Cube(domainSize); 
        Box patchBox = domainBox.coarsen(boxSize);
        std::vector<Point> patches;
        for (auto patch : patchBox)
        {
            if (patch != Point::Zeros()) { patches.push_back(patch); }
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
}
TEST(HDF5, ReadWritePatch)
{
    Box B = Box::Cube(8);
    BoxData<double, 3, HOST> data(B);
    forallInPlace_p(f_pointID, data);
    HDF5Handler h5;
    h5.writePatch(data, "TEST_WRITE_PATCH");
    //todo: implement readPatch
}
TEST(HDF5, ReadWriteLevel)
{
    int domainSize = 32;
    int boxSize = 8;
    int ghostSize = 3;

    double dx = 1.0/domainSize;
    Point offset(1,2,3,4,5,6);
    Point k(1,2,3,4,5,6);

    auto layout = testLayout(domainSize, Point::Ones(boxSize));
    layout.print();
    LevelBoxData<double, 3, HOST> writeData(layout, Point::Ones(ghostSize));
    LevelBoxData<double, 3, HOST> readData(layout, Point::Ones(ghostSize));
    writeData.initialize(f_phi, dx, k, offset);

    HDF5Handler h5;
    h5.writeLevel(writeData, "TEST_WRITE_LEVEL");
    h5.readLevel(readData, "TEST_WRITE_LEVEL");

    LevelBoxData<double, 3, HOST> errorData(layout, Point::Ones(ghostSize));
    errorData.setVal(0);
    double error = 0;
    for (auto index : layout)
    {
        auto& wi = writeData[index];
        auto& ri = readData[index];
        auto& ei = errorData[index];
        EXPECT_EQ(wi.box(), ri.box());

        wi.copyTo(ei);
        ei -= ri;
        error = max(error, ei.absMax());
    }
    EXPECT_LT(error, 1e-12);
}
TEST(HDF5, ReadWriteMBLevel)
{
    HDF5Handler h5;
    int domainSize = 4;
    int boxSize = 2;
    auto domain = buildXPoint(domainSize);
    Point boxSizeVect = Point::Ones(boxSize);
    std::vector<MBPoint> patches;
    std::vector<Point> boxSizes;
    for (BlockIndex bi = 0; bi < domain.numBlocks(); bi++)
    {
        boxSizes.push_back(Point::Ones(boxSize));
        for (Point pi : Box::Cube(domainSize / boxSize))
        {
            if (pi == Point::Zeros()) { continue; }
            patches.push_back(MBPoint(pi, bi));
        }
    }

    MBDisjointBoxLayout layout(domain, patches, boxSizes);
    MBLevelBoxData<double, DIM, HOST> data(layout, Point::Ones());
    data.initialize(f_MBPointID);
    h5.writeMBLevel(data, "TEST_WRITE_MB_LEVEL");
    
    MBDisjointBoxLayout newLayout;
    h5.readMBLayout(newLayout, domain.graphPtr(), "TEST_WRITE_MB_LEVEL");
    EXPECT_TRUE(newLayout.compatible(layout));
    for (auto iter : layout)
    {
        EXPECT_EQ(layout[iter], newLayout[iter]);
    }

    MBLevelBoxData<double, DIM, HOST> newData(newLayout, Point::Ones());
    newData.setVal(7);
    h5.readMBLevel(newData, "TEST_WRITE_MB_LEVEL");

    for (auto iter : layout)
    {
        auto& old_i = data[iter];
        auto& new_i = newData[iter];

        EXPECT_EQ(old_i.box(), new_i.box());
        BoxData<double, DIM, HOST> error(new_i.box());
        new_i.copyTo(error);
        error -= old_i;
        EXPECT_NEAR(error.absMax(), 0, 1e-12);
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
