#include <gtest/gtest.h>
#include "Proto.H"
#include "Lambdas.H"
#include "BoxOp_MBLaplace.H"
#include "MBMap_XPointRigid.H"

using namespace Proto;


class Data {
    public:
        Data(const MBDisjointBoxLayout& a_layout, int a_refRatio)
        {
            std::vector<Point> refRatios(a_layout.numBlocks(), Point::Ones(a_refRatio));
            Array<Point, DIM+1> ghost;
            ghost.fill(Point::Zeros());
            auto crseLayout = a_layout.coarsen(refRatios);
            //m_data.define(crseLayout, ghost);
            auto domain = a_layout.domain();
            auto boxSize = a_layout.boxSizes()[0];
            MBDisjointBoxLayout layout(domain, boxSize);
            m_data.define(layout, ghost);
        }
    
        MBLevelBoxData<double, 1, HOST>& get() { return m_data; }
    private:
        MBLevelBoxData<double, 1, HOST> m_data;
};
#if 0
TEST(MBMultigridTests, Temp) {
    int domainSize = 32;
    int boxSize = 16;
    int numBlocks = 5;
    int numLevels = 3;
    Point refRatio = Point::Ones(2);
    std::vector<Point> refRatios;
    for (int bi = 0; bi < numBlocks; bi++)
    {
        refRatios.push_back(refRatio);
    }
    auto domain = buildXPoint(domainSize, numBlocks);
    MBDisjointBoxLayout layout(domain, Point::Ones(boxSize));
    
    Data D(layout, 2);
    std::cout << &(D.get()) << std::endl;
    std::cout << &(D.get().layout()) << std::endl;
    std::cout << &(D.get().layout().partition()) << std::endl;
    std::cout << D.get().layout().numBlocks() << std::endl;
}
#endif

#if 1
TEST(MBMultigridTests, Construction) {
    int domainSize = 32;
    int boxSize = 16;
    int numBlocks = 5;
    int numLevels = 3;
    Point refRatio = Point::Ones(2);
    std::vector<Point> refRatios;
    for (int bi = 0; bi < numBlocks; bi++)
    {
        refRatios.push_back(refRatio);
    }
    auto domain = buildXPoint(domainSize, numBlocks);
    MBDisjointBoxLayout layout(domain, Point::Ones(boxSize));

    MBMultigrid<BoxOp_MBLaplace, MBMap_XPointRigid, double> mg(layout, refRatio, numLevels); 
    EXPECT_TRUE(mg.validate(layout, refRatios, numLevels));
}
#endif

int main(int argc, char *argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
#ifdef PR_MPI
    MPI_Init(&argc, &argv);
#endif
    PR_TIMER_SETFILE("MBMultigridTests_DIM" + to_string(DIM) + "_NProc" + to_string(numProc())
            + ".time.table");
    int result = RUN_ALL_TESTS();
    PR_TIMER_REPORT();
#ifdef PR_MPI
    MPI_Finalize();
#endif
    return result;
}
