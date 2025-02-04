#include <gtest/gtest.h>
#include "Proto.H"
#include "TestFunctions.H"
#include "BoxOp_MBLaplace.H"
#include "MBMap_XPointRigid.H"

using namespace Proto;

namespace {
    template< typename T, unsigned int C, MemType MEM, Centering CTR>
    void initializeData(MBAMRData<T,C,MEM,CTR>& data)
    {
        Array<double, DIM> offset{0,0,0,0,0,0};

        auto C2C = Stencil<double>::CornersToCells(4);
        for (int li = 0; li < numLevels; li++)
        {
            auto& levelLayout = layout[li];
            for (auto iter : levelLayout)
            {
                auto block = levelLayout.block(iter);
                auto& data_i = data[li][iter];
                Box b_i = C2C.domain(levelLayout[iter]).grow(data.ghost()[0]);
                BoxData<double, DIM> x_i(b_i.grow(PR_NODE));
                BoxData<double, 1> J_i(x_i.box());
                // compute coordinate / jacobian values from map
                map[li].apply(x_i, J_i, block);
                BoxData<double, 1> nodeData = forall<double, 1>(f_bell, x_i, offset);

                data_i |= C2C(nodeData);
            }
        }
    }
}

#if 1
TEST(MBAMR, AverageDown) {
    HDF5Handler h5;
    int domainSize = 16;
    int boxSize = 16;
    int numBlocks = 5;
    int numLevels = 2;
    int refRatio = 2;
    int ghostWidth = 1;
    Point ghostWidths = Point::Ones(ghostWidth);
  
    int N = 1;
    double error[N];
    for (int nn = 0; nn < N; nn++)
    {
        MBAMRLayout layout = telescopingXPointGrid(domainSize, numLevels, refRatio, boxSize, numBlocks);
        MBAMR<BoxOp_MBLaplace, MBMap_XPointRigid, double> amr(layout, Point::Ones(refRatio)); 
        MBAMRMap<MBMap_XPointRigid, HOST> map(layout, ghostWidths);
        for (int li = 0; li < numLevels; li++)
        {
            for (BlockIndex bi = 0; bi < numBlocks; bi++)
            {
                map[li][bi].setNumBlocks(numBlocks);
            }
        }
        MBAMRData<double, 1, HOST> phi(layout, ghostWidths);
        MBAMRData<double, 1, HOST> err(layout, ghostWidths);
        MBAMRData<double, 1, HOST> phi0(layout, ghostWidths);
        initializeData(phi);
        initializeData(phi0);
        err.setVal(0);
        
        auto& crse = phi[0]; 
        auto& fine = phi0[1];
        amr.averageDown(crse, fine, 1);
#if PR_VERBOSE > 0
        h5.writeMBAMRData({"phi"}, map, phi0, "MBAMRTests_AverageDown_PHI_INIT_%i", nn);
        h5.writeMBAMRData({"phi"}, map, phi, "MBAMRTests_AverageDown_PHI_AVG_%i", nn);
#endif
        auto& levelLayout = layout[0];
        error[nn] = 0;
        for (auto iter : levelLayout)
        {
            auto& phi0_i = phi0[0][iter];
            auto& phi_i = phi[0][iter];
            auto& err_i = err[0][iter];
            err_i.setVal(0);
            phi_i.copyTo(err_i);
            err_i -= phi0_i;
            error[nn] = max(error[nn], err_i.absMax());
        }
        Reduction<double, Max> rxn;
        rxn.reduce(&error[nn], 1);
        error[nn] = rxn.fetch();
#if PR_VERBOSE > 0
        std::cout << "Error (Max Norm): " << error[nn] << std::endl;
        h5.writeMBAMRData({"err"}, map, err, "MBAMRTests_AverageDown_ERR_%i", nn);
#endif
        domainSize *= 2;
    }

    for (int ii = 1; ii < N; ii++)
    {
        double rate = log(error[ii-1]/error[ii])/log(2.0);
        EXPECT_TRUE(error[ii] < 1e-10 || abs(4.0-rate) < 0.3);
#if PR_VERBOSE > 0
        std::cout << "Convergence Rate: " << rate << std::endl;
#endif
    }
}
#endif

#if 1
TEST(MBAMR, InterpBounds) {
    HDF5Handler h5;
    int domainSize = 16;
    int boxSize = 16;
    int numBlocks = 5;
    int numLevels = 2;
    int refRatio = 2;
    int numGhost = 1;
    Array<double, DIM> offset{0,0,0,0,0,0};
    Array<Point,DIM+1> ghost;
    ghost.fill(Point::Ones(numGhost));
  
    int N = 3;
    double error[N];
    for (int nn = 0; nn < N; nn++)
    {
        MBAMRLayout grid = telescopingXPointGrid(domainSize, numLevels, refRatio, boxSize, numBlocks);
        MBAMR<BoxOp_MBLaplace, MBMap_XPointRigid, double> amr(
                grid, Point::Ones(refRatio)); 

        MBAMRMap<MBMap_XPointRigid, HOST> map(grid, ghost);
        for (int li = 0; li < numLevels; li++)
        {
            for (BlockIndex bi = 0; bi < numBlocks; bi++)
            {
                map[li][bi].setNumBlocks(numBlocks);
            }
        }
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
        amr.interpBounds(crse, fine, 1);
#if PR_VERBOSE > 0
        h5.writeMBAMRData({"phi"}, map, phi0, "MBAMRTests_InterpBounds_PHI_INIT_%i", nn);
        h5.writeMBAMRData({"phi"}, map, phi, "MBAMRTests_InterpBounds_PHI_AVG_%i", nn);
#endif
        auto& layout = grid[0];
        error[nn] = 0;
        for (auto iter : layout)
        {
            auto& phi0_i = phi0[0][iter];
            auto& phi_i = phi[0][iter];
            auto& err_i = err[0][iter];
            err_i.setVal(0);
            phi_i.copyTo(err_i);
            err_i -= phi0_i;
            error[nn] = max(error[nn], err_i.absMax());
        }
        Reduction<double, Max> rxn;
        rxn.reduce(&error[nn], 1);
        error[nn] = rxn.fetch();
#if PR_VERBOSE > 0
        std::cout << "Error (Max Norm): " << error[nn] << std::endl;
        h5.writeMBAMRData({"err"}, map, err, "MBAMRTests_InterpBounds_ERR_%i", nn);
#endif
        domainSize *= 2;
    }

    for (int ii = 1; ii < N; ii++)
    {
        double rate = log(error[ii-1]/error[ii])/log(2.0);
        EXPECT_TRUE(error[ii] < 1e-10 || abs(4.0-rate) < 0.3);
#if PR_VERBOSE > 0
        std::cout << "Convergence Rate: " << rate << std::endl;
#endif
    }
}
#endif
#if 1
TEST(MBAMRTests, LaplaceXPoint) {
    HDF5Handler h5;
    
    // Constants
    typedef BoxOp_MBLaplace<double, MBMap_XPointRigid> OP;
    int domainSize = 16;
    int boxSize = 16;
    int numBlocks = 5;
    int numLevels = 1;
    Array<double, DIM> offset{0,0,0,0,0,0};
    int refRatio = 2;
   
    for (int nn = 0; nn < 1; nn++)
    {
        MBAMRLayout grid = telescopingXPointGrid(domainSize, numLevels, refRatio, boxSize);
        MBAMR<BoxOp_MBLaplace, MBMap_XPointRigid, double> amr(
            grid, Point::Ones(refRatio));

        MBAMRMap<MBMap_XPointRigid, HOST> map(grid, OP::ghost());
        for (int li = 0; li < numLevels; li++)
        {
            for (BlockIndex bi = 0; bi < numBlocks; bi++)
            {
                map[li][bi].setNumBlocks(numBlocks);
            }
        }
        MBAMRData<double, 1, HOST> phi(grid, OP::ghost());
        MBAMRData<double, 1, HOST> rhs(grid, Point::Zeros());

        auto C2C = Stencil<double>::CornersToCells(4);
        for (int li = 0; li < numLevels; li++)
        {
            auto& layout = grid[li];
            for (auto iter : layout)
            {
                auto block = layout.block(iter);
                auto& rhs_i = rhs[li][iter];
                Box b_i = C2C.domain(layout[iter]).grow(OP::ghost());
                BoxData<double, DIM> x_i(b_i.grow(PR_NODE));
                BoxData<double, 1> J_i(b_i);
                // compute coordinate / jacobian values from map
                map[li].apply(x_i, J_i, block);
                BoxData<double, 1> rhs_node = forall<double, 1>(f_bell, x_i, offset);

                rhs_i |= C2C(rhs_node);
            }
        }
        phi.setVal(0);
       
#if PR_VERBOSE > 0
        h5.writeMBAMRData({"rhs"}, map, rhs, "MBAMRTests_LaplaceXPoint_RHS_%i", nn);
        h5.writeMBAMRData({"phi"}, map, phi, "MBAMRTests_LaplaceXPoint_PHI_%i_0", nn);
#endif
        amr.solve(phi, rhs, 10, 1e-10);
#if PR_VERBOSE > 0
        h5.writeMBAMRData({"phi"}, map, phi, "MBAMRTests_LaplaceXPoint_PHI_%i_1", nn);
#endif

        domainSize *= 2;
        boxSize *= 2;
    }
}
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
