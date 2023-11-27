#include <gtest/gtest.h>
#include "Proto.H"
#include "Lambdas.H"
#include "BoxOp_MBLaplace.H"
#include "MBMap_XPointRigid.H"

using namespace Proto;

PROTO_KERNEL_START
void f_linearF(Var<double>& a_data, Var<double, DIM>& a_X, double a_slope, int a_comp)
{
    a_data(0) = a_X(a_comp)*a_slope;
}
PROTO_KERNEL_END(f_linearF, f_linear)

PROTO_KERNEL_START
void f_bellF(Var<double>& a_data, Var<double, DIM>& a_X, Array<double, DIM> a_X0)
{
    double p = 6.0;
    double s = 2.0;
    double r2 = 0;
    for (int ii = 0; ii < DIM; ii++)
    {
        double xi = (a_X(ii) - a_X0[ii])*s;
        r2 += xi*xi;
    }
    double r = sqrt(r2);
    a_data(0) = (r < M_PI/2.0) ? pow(cos(r),p) : 0.0;
}
PROTO_KERNEL_END(f_bellF, f_bell)

PROTO_KERNEL_START
void f_LbellF(Var<double>& a_data, Var<double, DIM>& a_X, Array<double, DIM> a_X0)
{
    double p = 6.0;
    double s = 2.0;
    double r2 = 0;
    for (int ii = 0; ii < DIM; ii++)
    {
        double xi = (a_X(ii) - a_X0[ii]);
        r2 += xi*xi;
    }
    double r = sqrt(r2);
    double sr = sin(s*r);
    double cr = cos(s*r);
    a_data(0) = 0;
    if (s*r < M_PI/2.0)
    {
        for (int dir = 0; dir < DIM; dir++)
        {
            double x = (a_X(dir) - a_X0[dir]);
            a_data(0) -= p*s*s*x*x*pow(cr,p)/(r*r);
            a_data(0) += p*(p-1)*s*s*x*x*pow(sr,2)*pow(cr,p-2)/(r*r);
            a_data(0) -= p*s*sr*pow(cr,p-1)/r;
            a_data(0) += p*s*x*x*sr*pow(cr,p-1)/(r*r*r);
            //    a_data(0) += p*(p-1)*pow(cr,p-2)*sr*sr*x*x/(r*r);
        //    a_data(0) += p*pow(cr,p-1)*sr*x*x/(r*r*r);
        //    a_data(0) -= p*pow(cr,p)*x*x/(r*r);
        //    a_data(0) -= p*pow(cr,p-1)*sr/r;
        }
    }
}
PROTO_KERNEL_END(f_LbellF, f_Lbell)

#if 1
TEST(MBMultigridTests, LaplaceXPoint) {
    HDF5Handler h5;

    typedef BoxOp_MBLaplace<double> OP;
    int domainSize = 32;
    int boxSize = 16;
    int numBlocks = XPOINT_NUM_BLOCKS;
    int numLevels = 1;
    double slope = 1.0;
    int comp = 0;
    Array<double, DIM> k{1,1,1,0,0,0};
    Array<double, DIM> offset{0,0,0,0,0,0};
    offset += 0.01;
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

    MBLevelBoxData<double, 1, HOST> phi(layout, OP::ghost());
    MBLevelBoxData<double, 1, HOST> lphi(layout, Point::Zeros());
    MBLevelBoxData<double, 1, HOST> rhs(layout, Point::Zeros());
    MBLevelBoxData<double, 1, HOST> res(layout, Point::Zeros());
    
    MBLevelMap<MBMap_XPointRigid, HOST> map(layout, OP::ghost());

    auto C2C = Stencil<double>::CornersToCells(4);
    for (auto iter : layout)
    {
        auto block = layout.block(iter);
        auto& phi_i = phi[iter];
        auto& rhs_i = rhs[iter];
        Box b_i = C2C.domain(layout[iter]).grow(OP::ghost());
        BoxData<double, DIM> x_i(b_i.grow(PR_NODE));
        BoxData<double, 1> J_i(b_i);
        map.apply(x_i, J_i, block);
#if 1
        /*
        BoxData<double, 1> phi0 = forall<double, 1>(f_linear, x_i, slope, comp);
        BoxData<double, 1> lphi(layout[iter]);
        lphi.setVal(0);
        */
        BoxData<double, 1> phi0 = forall<double, 1>(f_bell, x_i, offset);
        BoxData<double, 1> lphi = forall<double, 1>(f_Lbell, x_i, offset);

#else
        BoxData<double, 1> phi0 = forall_p<double, 1>(f_phiM, block, x_i, k, offset);
        BoxData<double, 1> lphi = forall_p<double, 1>(f_LphiM, block, x_i, k, offset);
#endif
        phi_i |= C2C(phi0);
        lphi *= J_i.absMax();
        rhs_i |= C2C(lphi);
    }
#if PR_VERBOSE > 0
    h5.writeMBLevel({"phi"}, map, phi, "PHI_0");
    h5.writeMBLevel({"rhs"}, map, rhs, "RHS_0");
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
    int result = RUN_ALL_TESTS();
#ifdef PR_MPI
    MPI_Finalize();
#endif
    return result;
}
