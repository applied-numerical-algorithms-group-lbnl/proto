#include <gtest/gtest.h>
#include "Proto.H"
#include "Lambdas.H"
#include "BoxOp_MBLaplace.H"
#include "MBMap_XPointRigid.H"

using namespace Proto;

#if 1
TEST(MBMultigridTests, LaplaceXPoint) {
    HDF5Handler h5;

    typedef BoxOp_MBLaplace<double> OP;
    int domainSize = 32;
    int boxSize = 16;
    int numBlocks = 5;
    int numLevels = 1;
    Array<double, DIM> k{1,1,1,0,0,0};
    Array<double, DIM> offset{0,0,0,0,0,0};
    offset += 0.1;
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

    Array<Point,DIM+1> ghost;
    Array<Point,DIM+1> noGhost;
    ghost.fill(OP::ghost());
    noGhost.fill(Point::Zeros());
    MBLevelBoxData<double, 1, HOST> phi(layout, ghost);
    MBLevelBoxData<double, 1, HOST> lphi(layout, noGhost);
    MBLevelBoxData<double, 1, HOST> rhs(layout, noGhost);
    MBLevelBoxData<double, 1, HOST> res(layout, noGhost);
    
    MBLevelMap<MBMap_XPointRigid, HOST> map(layout, ghost);

    auto C2C = Stencil<double>::CornersToCells(4);
    for (auto iter : layout)
    {
        auto block = layout.block(iter);
        auto& phi_i = phi[iter];
        Box b_i = C2C.domain(layout[iter]).grow(ghost[0]);
        BoxData<double, DIM> x_i(b_i);
        BoxData<double, 1> J_i(b_i);
        map.apply(x_i, J_i, block);
        BoxData<double, 1> phi0 = forall_p<double, 1>(f_phiM, block, x_i, k, offset);
        phi_i |= C2C(phi0);

        auto& rhs_i = rhs[iter];
        BoxData<double, 1> lphi = forall_p<double, 1>(f_LphiM, block, x_i, k, offset);
        lphi *= J_i.absMax(); //J is a constant
        //BoxData<double, 1> lphi_c = C2C(lphi, 1);
        //Operator::cellProduct(rhs_i, lphi_c, J_i);
        rhs_i |= C2C(lphi);
    }
    mg.op(numLevels-1)(lphi, phi);
    //mg.residual(res, phi, rhs, numLevels-1);
    

#if PR_VERBOSE > 0
    h5.writeMBLevel({"phi"}, map, phi, "PHI_");
#endif
    mg.solve(phi, rhs, 1, 1e-10);
#if PR_VERBOSE > 0
    h5.writeMBLevel({"phi"}, map, phi, "PHI_1");
#endif
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
