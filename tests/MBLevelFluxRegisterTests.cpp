#include <gtest/gtest.h>
#include "Proto.H"
#include "TestFunctions.H"
#include "MBMap_XPointRigid.H"

using namespace Proto;

// namespace {
//     <template template<MemType> typename MAP>
//     void flux_0(
//             FluxBoxData<double, 1>& a_flux,
//             MAP a_map,
//             MBIndex a_index)
//     {
//         auto& layout = a_map.layout();
//         Array<double, DIM> velocity = {1,1,1;};
//         double rho = 1.0;

//         auto block = layout.block(a_index);

//         // flux = rho*vd
//         for (int dd = 0; dd < DIM; dd++)
//         {
//             a_flux[dd].setVal(rho);
//             a_flux[dd] *= velocity[dd];
//         }
//     }
//}
#if PR_MMB
#if PR_AMR
TEST(MBLevelFluxRegister, Reflux_XPoint) {
    constexpr int NUM_COMPS = 2;
    int domainSize = 16;
    int boxSize = 16;
    Point boxSizes = Point::Ones(boxSize);
    int numBlocks = MB_MAP_XPOINT_NUM_BLOCKS;
    int numLevels = 2;
    int refRatio = 2;
    Point refRatios = Point::Ones(refRatio);
    int ghostWidth = 2;
    Point ghostWidths = Point::Ones(ghostWidth));
    
    auto grid = telescopingXPointGrid(domainSize, numLevels, refRatio, boxSize);
    MBAMRMap<MBMap_XPointRigid, HOST> map(grid, ghost);
    MBAMRData<double, 1, HOST> data(grid, ghost);

    Array<double, DIM> gridSpacing = Point::Ones();
    gridSpacing *= (1.0/domainSize);

    MBLevelFluxRegister<double, NUM_COMPS, HOST>(grid[0], grid[1], refRatios, gridSpacing);
    
}

#if DIM == 3
TEST(MBLevelFluxRegister, Reflux_CubedSphere) {
   
    int domainSize = 32;
    int boxSize = 16;
    int ghostSize = 2;

}
#endif

#endif
#endif
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
