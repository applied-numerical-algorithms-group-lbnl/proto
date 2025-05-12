#include <gtest/gtest.h>
#include "Proto.H"
#include "TestFunctions.H"
#include "BoxOp_MBLaplace.H"
#include "MBMap_XPointRigid.H"
#include "MBMap_Identity.H"

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
#if 0
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
#endif
#if 0
TEST(MBMultigridTests, RelaxSingleIter)
{
    #if PR_VERBOSE > 0
    HDF5Handler h5;
    #endif

    int domainSize = 16;
    int boxSize = 8;
    int numBlocks = 4;
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
        mg.map(lvl).initialize();
    }

    MBLevelBoxData<double, 1, HOST> phi(layout, OP::ghost());
    MBLevelBoxData<double, 1, HOST> err(layout, Point::Zeros());
    MBLevelBoxData<double, 1, HOST> rhs(layout, Point::Zeros());
    MBLevelBoxData<double, 1, HOST> res(layout, Point::Zeros());

    initializeData(phi, rhs, mg.map());
    phi.setVal(0);
    res.setVal(0);

    mg.residual(res, phi, rhs);
    mg.relax(phi, rhs, 1);

    err.setVal(0);
    for (auto index : layout)
    {
        // error is phi1 - (-lambda*res)
        res[index] *= mg.lambda();
        err[index] += res[index];
        err[index] += phi[index];
    }

    EXPECT_LT(err.absMax(), 1e-12);

#if PR_VERBOSE > 0
    h5.writeMBLevel({"phi"}, mg.map(), phi, "TEST_RELAX_PHI");
    h5.writeMBLevel({"lambda*res"}, mg.map(), res, "TEST_RELAX_LAMBDA_RES");
    h5.writeMBLevel({"err"}, mg.map(), err, "TEST_RELAX_ERR");
#endif

}
#endif
#if 1
TEST(MBMultigridTests, RelaxConvergence)
{
    #if PR_VERBOSE > 0
    HDF5Handler h5;
    #endif

    int domainSize = 16;
    int boxSize = 16;
    int numBlocks = 4;
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
        mg.map(lvl).initialize();
    }
    mg.interpOp().writeFootprint("MG_INTERP");
    MBLevelBoxData<double, 1, HOST> phi(layout, OP::ghost());
    MBLevelBoxData<double, 1, HOST> rhs(layout, Point::Zeros());
    MBLevelBoxData<double, 1, HOST> res(layout, Point::Zeros());

    initializeData(phi, rhs, mg.map());
    phi.setVal(0);
    res.setVal(0);

    mg.residual(res, phi, rhs);
    mg.relax(phi, rhs, 2);
}
#endif
#if 0
TEST(MBMultigridTests, AverageDown)
{
    #if PR_VERBOSE > 0
    HDF5Handler h5;
    #endif

    int domainSize = 64;
    int boxSize = 32;
    int numBlocks = 5;
    int numLevels = 2;
    Point refRatio = Point::Ones(2);

    int numIter = 2;
    double error[numIter];
    for (int nn = 0; nn < numIter; nn++)
    {

        auto domain = buildXPoint(domainSize, numBlocks);
        MBDisjointBoxLayout layout(domain, Point::Ones(boxSize));

        MBMultigrid<BoxOp_MBLaplace, MBMap_XPointRigid, double> mg(layout, refRatio, numLevels);
        for (int lvl = 0; lvl < numLevels; lvl++)
        {
            for (BlockIndex bi = 0; bi < numBlocks; bi++)
            {
                mg.map(lvl)[bi].setNumBlocks(numBlocks);
            }
            mg.map(lvl).initialize();
        }

        MBLevelBoxData<double, 1, HOST> phiAvg(mg.map(0).layout(), Point::Zeros());
        MBLevelBoxData<double, 1, HOST> phiCrse(mg.map(0).layout(), Point::Zeros());
        MBLevelBoxData<double, 1, HOST> phiErr(mg.map(0).layout(), Point::Zeros());
        MBLevelBoxData<double, 1, HOST> phiFine(mg.map(1).layout(), Point::Zeros());
        MBLevelBoxData<double, 1, HOST> rhsCrse(mg.map(0).layout(), Point::Zeros());
        MBLevelBoxData<double, 1, HOST> rhsFine(mg.map(1).layout(), Point::Zeros());

        initializeData(phiCrse, rhsCrse, mg.map(0));
        initializeData(phiFine, rhsFine, mg.map(1));

        phiAvg.setVal(0);
        phiErr.setVal(0);
        mg.averageDown(phiAvg, phiFine, 1);
        phiErr.increment(phiAvg);
        phiErr.increment(phiCrse, -1);

        error[nn] = phiErr.absMax();

#if PR_VERBOSE > 0
        std::cout << "error: " << phiErr.absMax() << std::endl;
        h5.writeMBLevel({"phi"}, mg.map(0), phiAvg, "TEST_AVGDOWN_PHI_AVG_%i",nn);
        h5.writeMBLevel({"phi"}, mg.map(1), phiFine, "TEST_AVGDOWN_PHI_FIN_%i",nn);
        h5.writeMBLevel({"phi"}, mg.map(0), phiCrse, "TEST_AVGDOWN_PHI_CRS_%i",nn);
        h5.writeMBLevel({"phi"}, mg.map(0), phiErr, "TEST_AVGDOWN_PHI_ERR_%i",nn);
#endif
        domainSize *= 2;
        boxSize *= 2;
    }
    for (int nn = 1; nn < numIter; nn++)
    {
        double rate = log(error[nn - 1] / error[nn]) / log(2.0);
        EXPECT_LT(abs(4-rate), 0.1);
#if PR_VERBOSE > 0
        std::cout << "rate: " << rate << std::endl;
#endif
    }
}
#endif
#if 0
TEST(MBMultigridTests, LaplaceIdentity) {
    #if PR_VERBOSE > 0
    HDF5Handler h5;
    #endif

    int domainSize = 32;
    int boxSize = 16;
    int numLevels = log(domainSize)/log(2.0);
    Point refRatio = Point::Ones(2);
   
    auto domain = buildIdentity(domainSize);
    MBDisjointBoxLayout layout(domain, Point::Ones(boxSize));

    MBMultigrid<BoxOp_MBLaplace, MBMap_Identity, double> mg(layout, refRatio, numLevels);


    MBLevelBoxData<double, 1, HOST> phi(layout, OP::ghost());
    MBLevelBoxData<double, 1, HOST> rhs(layout, Point::Zeros());
    MBLevelBoxData<double, 1, HOST> res(layout, Point::Zeros());

    MBLevelMap<MBMap_Identity, HOST> map(layout, OP::ghost());

    initializeData(phi, rhs, map);
    
    mg.op(numLevels-1)(rhs, phi);
    phi.setVal(0);
#if PR_VERBOSE > 0
    h5.writeMBLevel({"rhs"}, map, rhs, "RHS_LAPLACE_IDENTITY");
    h5.writeMBLevel({"phi"}, map, phi, "PHI_LAPLACE_IDENTITY");
    h5.writeMBLevel({"residual"}, map, res, "RES_LAPLACE_IDENTITY");
#endif
    mg.solve(phi, rhs, 10, 1e-10);
}
#endif
#if 0
TEST(MBMultigridTests, LaplaceXPoint) {
    #if PR_VERBOSE > 0
    HDF5Handler h5;
    #endif

    int domainSize = 16;
    int boxSize = 8;
    int numBlocks = 4;
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
    res.setVal(0);
#if PR_VERBOSE > 0
    h5.writeMBLevel({"rhs"}, map, rhs, "RHS_LAPLACE_XPOINT");
    h5.writeMBLevel({"phi"}, map, phi, "PHI_LAPLACE_XPOINT");
    h5.writeMBLevel({"residual"}, map, res, "RES_LAPLACE_XPOINT");
#endif
    mg.solve(phi, rhs, 10, 1e-10);
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
