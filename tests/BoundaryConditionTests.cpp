#include <gtest/gtest.h>
#include "Proto.H"
#include "TestFunctions.H"

using namespace Proto;

#if DIM==2
TEST(BoundaryCondition, ConstDirichlet)
{
    int domainSize = 8;
    Array<double, DIM> dx;
    dx.fill(1.0/domainSize);
    Box B0 = Box::Cube(domainSize);
    Box BX = B0.extrude(Point::X());
    Box BY = B0.extrude(Point::Y());

    Face fx_lo = Face(0,Side::Lo);
    Face fx_hi = Face(0,Side::Hi);
    Face fy_lo = Face(1,Side::Lo);
    Face fy_hi = Face(1,Side::Hi);

    Array<BoxData<double>,DIM> flux;
    flux[0].define(BX);
    flux[1].define(BY);
    BoxData<double> phi(B0);

#if 0
    phi.setVal(0.0);
#else
    auto C2C = Stencil<double>::CornersToCells(4);
    Box BXY = C2C.domain(B0);
    auto phi0 = forall_p<double>(
            [] PROTO_LAMBDA (Point p, Var<double>& v, Array<double,DIM> dx)
            {
                Array<double, DIM> x;
                v(0) = 1.0;
                for (int ii = 0; ii < DIM; ii++)
                {
                    x[ii] = p[ii]*dx[ii];
                    v(0) *= cos(2.0*M_PI*x[ii]);
                }
            }, BXY, dx);
    phi |= C2C(phi0);
#endif
    for (int dd = 0; dd < DIM; dd++) { flux[dd].setVal(0); }

    pr_out() << "Phi: " << std::endl;
    phi.printData();
    pr_out() << "FX: " << std::endl;
    flux[0].printData();
    pr_out() << "FY: " << std::endl;
    flux[1].printData();

    BoundaryCondition::Dirichlet(flux, phi, 1.0, fx_lo, B0, dx);
    pr_out() << "FX (after lo BC): " << std::endl;
    flux[0].printData();
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
