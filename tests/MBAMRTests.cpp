#include <gtest/gtest.h>
#include "Proto.H"
#include "Lambdas.H"
#include "BoxOp_MBLaplace.H"
#include "MBMap_XPointRigid.H"

using namespace Proto;

TEST(MBAMR, Construction) {
    HDF5Handler h5;
    int domainSize = 16;
    int boxSize = 16;
    int numBlocks = XPOINT_NUM_BLOCKS;
    int numLevels = 2;
    int refRatio = 2;
    int numGhost = 1;
    Array<double, DIM> offset{0,0,0,0,0,0};
    Array<Point,DIM+1> ghost;
    ghost.fill(Point::Ones(numGhost));
    
    /*
    auto domain = buildXPoint(domainSize);
    std::vector<Point> boxSizeVect(numBlocks, Point::Ones(boxSize));
    std::vector<Point> refRatios(numLevels-1, Point::Ones(refRatio));

    auto coarsePatches = domain.patches(Point::Ones(boxSize));
    MBAMRGrid grid(domain, coarsePatches, boxSizeVect, refRatios);
    for (int li = 1; li < numLevels; li++)
    {
        auto domain = grid[li].domain();
        std::vector<MBPatchID_t> patches;
        Box patchDomain = Box::Cube(domainSize).refine(pow(refRatio,li)).coarsen(boxSize);
        Point patch = patchDomain.high();
        for (int bi = 0; bi < domain.size(); bi++)
        {
            patches.push_back(MBPatchID_t(patch, bi));
        }
        grid[li].define(domain, patches, boxSizeVect);
    }
    */

    MBAMRGrid grid = telescopingXPointGrid(domainSize, numLevels, refRatio, boxSize);
    for (int li = 0; li < numLevels; li++)
    {
        for (auto iter : grid[li])
        {
            std::cout << "p: " << grid[li].point(iter) << " | b: " << grid[li].block(iter) << " | l: " << li << std::endl;
        }
    }


    MBAMR<BoxOp_MBLaplace, MBMap_XPointRigid, double> amr(
            grid, Point::Ones(refRatio)); 

    MBAMRMap<MBMap_XPointRigid, HOST> map(grid, ghost);
    MBAMRData<double, 1, HOST> phi(grid, ghost);
    MBAMRData<double, 1, HOST> err(grid, ghost);
    MBAMRData<double, 1, HOST> phi0(grid, ghost);
    
    auto C2C = Stencil<double>::CornersToCells(4);
    for (int li = 0; li < numLevels; li++)
    {
        auto& layout = grid[li];
        for (auto iter : layout)
        {
            auto block = layout.block(iter);
            auto& phi0_i = phi0[li][iter];
            auto& phi_i = phi[li][iter];
            Box b_i = C2C.domain(layout[iter]).grow(numGhost);
            BoxData<double, DIM> x_i(b_i.grow(PR_NODE));
            BoxData<double, 1> J_i(b_i);
            // compute coordinate / jacobian values from map
            map[li].apply(x_i, J_i, block);
            BoxData<double, 1> phi_node = forall<double, 1>(f_bell, x_i, offset);
            
            phi_i |= C2C(phi_node);
            phi0_i |= C2C(phi_node);
        }

    }
    err.setVal(0);
    auto& crse = phi[0]; 
    auto& fine = phi0[1];
    amr.averageDown(crse, fine, 1);
#if PR_VERBOSE > 0
    h5.writeMBAMRData({"phi"}, map, phi0, "PHI_0");
    h5.writeMBAMRData({"phi"}, map, phi, "PHI_1");
#endif
    auto& layout = grid[0];
    for (auto iter : layout)
    {
        auto& phi0_i = phi0[0][iter];
        auto& phi_i = phi[0][iter];
        auto& err_i = err[0][iter];
        err_i.setVal(0);
        phi_i.copyTo(err_i);
        err_i -= phi0_i;
    }
#if PR_VERBOSE > 0
    h5.writeMBAMRData({"err"}, map, err, "ERR");
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
