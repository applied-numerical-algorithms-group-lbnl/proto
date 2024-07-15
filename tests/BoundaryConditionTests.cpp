#include <gtest/gtest.h>
#include "Proto.H"
#include "Lambdas.H"

using namespace Proto;

TEST(BoundaryCondition, ConstDirichlet)
{
    int domainSize = 8;
    double dx = 1.0/domainSize;
    Box B0 = Box::Cube(domainSize);
    Box BX = B0.extrude(Point::X());
    Box BY = B0.extrude(Point::Y());

    Face fx_lo = Face(0,Side::Lo);
    Face fx_hi = Face(0,Side::Hi);
    Face fy_lo = Face(1,Side::Lo);
    Face fy_hi = Face(1,Side::Hi);

    BoxData<double> flux_x(BX);
    BoxData<double> flux_y(BY);
    BoxData<double> phi(B0);

#if 0
    phi.setVal(0.0);
#else
    auto C2C = Stencil<double>::CornersToCells(4);
    Box BXY = C2C.domain(B0);
    auto phi0 = forall_p<double>(
            [] PROTO_LAMBDA (Point p, Var<double>& v, double dx)
            {
                Array<double, DIM> x;
                v(0) = 1.0;
                for (int ii = 0; ii < DIM; ii++)
                {
                    x[ii] = p[ii]*dx;
                    v(0) *= cos(2.0*M_PI*x[ii]);
                }
            }, BXY, dx);
    phi |= C2C(phi0);
#endif
    flux_x.setVal(0);
    flux_y.setVal(0);

    pr_out() << "Phi: " << std::endl;
    phi.printData();
    pr_out() << "FX: " << std::endl;
    flux_x.printData();
    pr_out() << "FY: " << std::endl;
    flux_y.printData();

    BoundaryCondition<double>::Dirichlet(flux_x, phi, B0, 1.0, dx, fx_lo);
    pr_out() << "FX (after lo BC): " << std::endl;
    flux_x.printData();
    /*
    BoundaryCondition<double>::Dirichlet(flux_x, phi, B0, dx, 1.0, fx_hi);
    pr_out() << "FX (after hi BC): " << std::endl;
    flux_x.printData();
    BoundaryCondition<double>::Dirichlet(flux_y, phi, B0, dx, 1.0, fy_lo);
    pr_out() << "FY (after lo BC): " << std::endl;
    flux_y.printData();
    BoundaryCondition<double>::Dirichlet(flux_y, phi, B0, dx, 1.0, fy_hi);
    pr_out() << "FY (after hi BC): " << std::endl;
    flux_y.printData();
*/
    BC_ConstDirichlet<double> bc_x_lo(fx_lo, 1.0, dx);
  //  BC_ConstDirichlet<double> bc_x_hi(fx_hi, 1.0, dx);
   // BC_ConstDirichlet<double> bc_y_lo(fy_lo, 1.0, dx);
   // BC_ConstDirichlet<double> bc_y_hi(fy_hi, 1.0, dx);

    flux_x.setVal(0);
    bc_x_lo.apply(flux_x, phi, B0);
    pr_out() << "FX (after lo BC): " << std::endl;
    flux_x.printData();
}

#if 0
TEST(BoundaryCondition, ConstDirichlet)
{
    int domainSize = 8;
    double dx = 1.0/domainSize;
    Box B0 = Box::Cube(domainSize);
    Box BX = B0.extrude(Point::X());
    Box BY = B0.extrude(Point::Y());
    auto C2C = Stencil<double>::CornersToCells(4);
    Box BXY = C2C.domain(B0);
    auto phi0 = forall_p<double>(
            [] PROTO_LAMBDA (Point p, Var<double>& v, double dx)
            {
                Array<double, DIM> x;
                v(0) = 1.0;
                for (int ii = 0; ii < DIM; ii++)
                {
                    x[ii] = p[ii]*dx;
                    v(0) *= cos(2.0*M_PI*x[ii]);
                }
            }, BXY, dx);
    BoxData<double> phi = C2C(phi0);
    phi0.printData();
    phi.printData();
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
