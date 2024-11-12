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
            patches.push_back(patch);
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
TEST(HDF5, ReadMBLevel)
{
    HDF5Handler h5;
    int domainSize = 4;
    int boxSize = 2;
    auto domain = buildXPoint(domainSize);
    Point boxSizeVect = Point::Ones(boxSize);
    std::vector<MBPatchID_t> patches;
    std::vector<Point> boxSizes;
    for (BlockIndex bi = 0; bi < domain.numBlocks(); bi++)
    {
        boxSizes.push_back(Point::Ones(boxSize));
        for (Point pi : Box::Cube(domainSize / boxSize))
        {
            if (pi == Point::Zeros()) { continue; }
            patches.push_back(std::make_pair(pi, bi));
        }
    }

    MBDisjointBoxLayout layout(domain, patches, boxSizes);
    MBLevelBoxData<double, DIM, HOST> data(layout, Point::Ones());
    data.initialize(f_MBPointID);
    h5.writeMBLevel(data, "DATA_0");
    
    MBDisjointBoxLayout newLayout;
    h5.readMBLayout(newLayout, domain.graphPtr(), "DATA_0");
    EXPECT_TRUE(newLayout.compatible(layout));
    for (auto iter : layout)
    {
        EXPECT_EQ(layout[iter], newLayout[iter]);
    }

    MBLevelBoxData<double, DIM, HOST> newData(newLayout, Point::Ones());
    newData.setVal(7);
    h5.readMBLevel(newData, "DATA_0");

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
