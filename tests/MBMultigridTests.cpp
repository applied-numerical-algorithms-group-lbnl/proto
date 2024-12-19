#include <gtest/gtest.h>
#include "Proto.H"
#include "TestFunctions.H"
#include "BoxOp_MBLaplace.H"
#include "MBMap_XPointRigid.H"

using namespace Proto;

PROTO_KERNEL_START
void f_linearF(Var<double>& a_data, Var<double, DIM>& a_X, double a_slope, int a_comp)
{
    a_data(0) = a_X(a_comp)*a_slope;
}
PROTO_KERNEL_END(f_linearF, f_linear)

#if 1
TEST(MBMultigridTests, LaplaceXPoint) {
    HDF5Handler h5;
    
    // Constants
    typedef BoxOp_MBLaplace<double, MBMap_XPointRigid> OP;
    int domainSize = 16;
    int boxSize = 8;
    int numBlocks = MB_MAP_XPOINT_NUM_BLOCKS;
    int numLevels = 1;
    double slope = 1.0;
    int comp = 0;
    Array<double, DIM> exp{4,4,0,0,0,0};
    Array<double, DIM> offset{0,0,0,0,0,0};
    Point refRatio = Point::Ones(2);
    std::vector<Point> refRatios;
    for (int bi = 0; bi < numBlocks; bi++)
    {
        refRatios.push_back(refRatio);
    }
   
    for (int nn = 0; nn < 1; nn++)
    {
        // Create Layout
        auto domain = buildXPoint(domainSize);
        MBDisjointBoxLayout layout(domain, Point::Ones(boxSize));

        // Create Multigrid Operator
        MBMultigrid<BoxOp_MBLaplace, MBMap_XPointRigid, double> mg(layout, refRatio, numLevels); 
        EXPECT_TRUE(mg.validate(layout, refRatios, numLevels));

#if PR_VERBOSE > 0
        h5.writeMBLevel({"x", "y"}, mg.map(numLevels-1), mg.map(numLevels-1).map(), "MAP_N%i",nn);
#endif        

        // Declare Data Holders
        MBLevelBoxData<double, 1, HOST> phi(layout, OP::ghost());
        MBLevelBoxData<double, 1, HOST> rhs(layout, Point::Zeros());
        MBLevelBoxData<double, 1, HOST> res(layout, Point::Zeros());

        // Build the (finest level) map
        MBLevelMap<MBMap_XPointRigid, HOST> map(layout, OP::ghost());
        
        // Initialize Data
        auto C2C = Stencil<double>::CornersToCells(4);
        double J0 = 0;
        for (auto iter : layout)
        {
            auto block = layout.block(iter);
            auto& phi_i = phi[iter];
            auto& rhs_i = rhs[iter];
            Box b_i = C2C.domain(layout[iter]).grow(OP::ghost());
            BoxData<double, DIM> x_i(b_i.grow(PR_NODE));
            BoxData<double, 1> J_i(b_i);
            // compute coordinate / jacobian values from map
            map.apply(x_i, J_i, block);
            BoxData<double, 1> phi_node = forall<double, 1>(f_bell, x_i, offset);
            BoxData<double, 1> lphi_node = forall<double, 1>(f_Lbell, x_i, offset);
            //BoxData<double, 1> phi_node = forall<double, 1>(f_poly, x_i, exp, offset);
            //BoxData<double, 1> lphi_node = forall<double, 1>(f_Lpoly, x_i, exp, offset);
            
            phi_i |= C2C(phi_node);
            rhs_i |= C2C(lphi_node);
            J0 = J_i.absMax();
        }

        mg.op(numLevels-1)(rhs, phi);
        phi.setVal(0);
#if PR_VERBOSE > 0
        h5.writeMBLevel({"rhs"}, map, rhs, "RHS_0_N%i",nn);
        h5.writeMBLevel({"phi"}, map, phi, "PHI_0_N%i",nn);
        h5.writeMBLevel({"residual"}, map, res, "RES_0_N%i",nn);
#endif
        mg.solve(phi, rhs, 20, 1e-10);
#if PR_VERBOSE > 0
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
