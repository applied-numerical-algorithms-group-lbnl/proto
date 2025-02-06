#include <gtest/gtest.h>
#include "Proto.H"
#include "TestFunctions.H"
#include "BoxOp_MBLaplace.H"
#include "MBMap_XPointRigid.H"

using namespace Proto;

namespace
{
    typedef BoxOp_MBLaplace<double, MBMap_XPointRigid> OP;

    PROTO_KERNEL_START
    void f_linearF(Var<double> &a_data, Var<double, DIM> &a_X, double a_slope, int a_comp)
    {
        a_data(0) = a_X(a_comp) * a_slope;
    }
    PROTO_KERNEL_END(f_linearF, f_linear)

    template <typename T, unsigned int C, MemType MEM, Centering CTR, template<MemType> typename MAP>
    void initializeData(
        MBLevelBoxData<T, C, MEM, CTR> &phi,
        MBLevelBoxData<T, C, MEM, CTR> &lphi,
        MBLevelMap<MAP, MEM>& map)
    {
        auto& layout = phi.layout();
        double dx = 0;
        Array<double, DIM> offset{dx,dx,0,0,0,0};
        auto C2C = Stencil<double>::CornersToCells(4);
        for (auto iter : layout)
        {
            auto block = layout.block(iter);
            auto &phi_i = phi[iter];
            auto &rhs_i = lphi[iter];
            Box b_i = C2C.domain(layout[iter]).grow(OP::ghost());
            BoxData<double, DIM> x_i(b_i.grow(PR_NODE));
            BoxData<double, 1> J_i(x_i.box());
            // compute coordinate / jacobian values from map
            map.apply(x_i, J_i, block);
            BoxData<double, 1> phi_node = forall<double, 1>(f_bell, x_i, offset);
            BoxData<double, 1> lphi_node = forall<double, 1>(f_Lbell, x_i, offset);

            phi_i |= C2C(phi_node);
            rhs_i |= C2C(lphi_node);
        }
    }
}

TEST(MBMultigridTests, Construction) 
{
    int domainSize = 16;
    int boxSize = 8;
    int numBlocks = 5;
    int numLevels = log(domainSize)/log(2.0);
    Point refRatio = Point::Ones(2);

    auto domain = buildXPoint(domainSize, numBlocks);
    MBDisjointBoxLayout layout(domain, Point::Ones(boxSize));
    MBMultigrid<BoxOp_MBLaplace, MBMap_XPointRigid, double> mg(layout, refRatio, numLevels);

    std::vector<Point> refRatios;
    for (int bi = 0; bi < numBlocks; bi++)
    {
        refRatios.push_back(refRatio);
    }
    EXPECT_TRUE(mg.validate(layout, refRatios, numLevels));
}

TEST(MBMultigridTests, Residual)
{
    #if PR_VERBOSE > 0
    HDF5Handler h5;
    #endif

    int domainSize = 16;
    int boxSize = 8;
    int numBlocks = 5;
    int numLevels = 1;
    Point refRatio = Point::Ones(2);

    auto domain = buildXPoint(domainSize, numBlocks);
    MBDisjointBoxLayout layout(domain, Point::Ones(boxSize));

    MBMultigrid<BoxOp_MBLaplace, MBMap_XPointRigid, double> mg(layout, refRatio, numLevels);
    for (int lvl = 0; lvl < numLevels; lvl++)
    {
        for (BlockIndex bi = 0; bi < numBlocks; bi++)
        {
            mg.map(lvl)[bi].setNumBlocks(numBlocks);
        }
    }

    MBLevelBoxData<double, 1, HOST> phi(layout, OP::ghost());
    MBLevelBoxData<double, 1, HOST> lphi(layout, Point::Zeros());
    MBLevelBoxData<double, 1, HOST> rhs(layout, Point::Zeros());
    MBLevelBoxData<double, 1, HOST> res(layout, Point::Zeros());
    MBLevelBoxData<double, 1, HOST> err(layout, Point::Zeros());

    initializeData(phi, rhs, mg.map());

    lphi.setVal(0);
    phi.exchange();
    mg.interpOp().apply(phi, phi);
    mg.op()(lphi, phi);

    mg.residual(res, phi, rhs);

    err.setVal(0);
    for (auto index : layout)
    {
        err[index] += res[index];
        err[index] -= rhs[index];
        err[index] += lphi[index];
    }

    EXPECT_LT(err.absMax(), 1e-12);
#if PR_VERBOSE > 0
    h5.writeMBLevel({"phi"}, mg.map(), phi, "TEST_RES_PHI");
    h5.writeMBLevel({"lphi"}, mg.map(), lphi, "TEST_RES_LPHI");
    h5.writeMBLevel({"rhs"}, mg.map(), rhs, "TEST_RES_RHS");
    h5.writeMBLevel({"residual"}, mg.map(), res, "TEST_RES_RES");
    h5.writeMBLevel({"error"}, mg.map(), err, "TEST_RES_ERR");
#endif
}

#if 1
TEST(MBMultigridTests, LaplaceXPoint) {
    #if PR_VERBOSE > 0
    HDF5Handler h5;
    #endif

    int domainSize = 16;
    int boxSize = 8;
    int numBlocks = 5;
    int numLevels = log(domainSize)/log(2.0);
    Point refRatio = Point::Ones(2);
   
    auto domain = buildXPoint(domainSize, numBlocks);
    MBDisjointBoxLayout layout(domain, Point::Ones(boxSize));

    MBMultigrid<BoxOp_MBLaplace, MBMap_XPointRigid, double> mg(layout, refRatio, numLevels);
    for (int lvl = 0; lvl < numLevels; lvl++)
    {
        for (BlockIndex bi = 0; bi < numBlocks; bi++)
        {
            mg.map(lvl)[bi].setNumBlocks(numBlocks);
        }
    }

    MBLevelBoxData<double, 1, HOST> phi(layout, OP::ghost());
    MBLevelBoxData<double, 1, HOST> rhs(layout, Point::Zeros());
    MBLevelBoxData<double, 1, HOST> res(layout, Point::Zeros());

    MBLevelMap<MBMap_XPointRigid, HOST> map(layout, OP::ghost());
    for (BlockIndex bi = 0; bi < numBlocks; bi++)
    {
        map[bi].setNumBlocks(numBlocks);
    }

    initializeData(phi, rhs, map);
    
    mg.op(numLevels-1)(rhs, phi);
    phi.setVal(0);
#if PR_VERBOSE > 0
    h5.writeMBLevel({"rhs"}, map, rhs, "RHS_0");
    h5.writeMBLevel({"phi"}, map, phi, "PHI_0");
    h5.writeMBLevel({"residual"}, map, res, "RES_0");
#endif
    mg.solve(phi, rhs, 20, 1e-10);

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
