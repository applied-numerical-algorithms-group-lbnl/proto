#include <gtest/gtest.h>
#include "Proto.H"
#include "TestFunctions.H"

using namespace Proto;



#if DIM==2
TEST(BoundaryCondition, ConstDirichletFlux)
{
    #if 1
    writeColorText(RED, "TEST DISABLED | This test should be either fixed or deleted");
    #else
    #if PR_VERBOSE > 0
    HDF5Handler h5;
    #endif
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
    #endif
}
#endif
namespace {
    PROTO_KERNEL_START
    void _f_initInterior(Point& pi, Var<double, 1, MEMTYPE_DEFAULT>& data, BlockIndex block, Array<double, DIM> dx)
    {
        Array<double, DIM> x(pi);
        x += 0.5;
        x *= dx;
        data(0) = 1.0;
        for (int ii = 0; ii < DIM; ii++)
        {
            data(0) *= sin(x[ii]); 
        }
    }
    PROTO_KERNEL_END(_f_initInterior, f_initInterior);
}

TEST(BoundaryCondition, ConstDirichletGhost)
{
    #if 1
    writeColorText(RED, "TEST DISABLED | This test sometimes fails when building with RELEASE. Should be investigated.");
    #else
    #if PR_VERBOSE > 0
    HDF5Handler h5;
    #endif
    int domainSize = 16;
    int boxSize = 8;
    int numIter = 3;
    int ghostSize = 2;
    constexpr int numBlocks = 5;

    double err[numIter];
    for (int nn = 0; nn < numIter; nn++)
    {
        auto domain = XPoint<numBlocks>::Domain(domainSize, boxSize);
        MBDisjointBoxLayout layout(domain, Point::Ones(boxSize));
        MBLevelMap<MBMap_XPointRigid<numBlocks, HOST>, HOST> map;
        map.define(layout, Point::Ones(ghostSize));
        MBInterpOp interp(map);
        MBLevelBoxData<double, 1, HOST> data(layout, Point::Ones(ghostSize));

        ConvolveOp<double> C;
        Array<double, DIM> dx;
        dx.fill(0.5*M_PI/domainSize);
        for (auto iter : layout)
        {
            auto [B0, B1] = C.domains(layout[iter].grow(ghostSize));
            auto tmp0 = forall_p<double, 1>(f_initInterior, B1, layout.block(iter), dx);
            C(data[iter], tmp0, tmp0);
        }
        data.setBoundary(1.0);

        #if PR_VERBOSE > 0
        HDF5Handler h5;
        h5.writeMBLevel(data, "DIRICHLET_GHOST_DATA_0");
        #endif

        Box domainBox = Box::Cube(domainSize);
        Box xBoundBox = domainBox.adjacent(-Point::X(),1);
        Box yBoundBox = domainBox.adjacent(-Point::Y(),1);
        err[nn] = 0;
        for (auto iter : layout)
        {
            auto& data_i = data[iter];
            Face fx(0, Side::Lo);
            Face fy(1,Side::Lo);
            BoundaryCondition::DirichletFillGhost<double>(data_i, 0, fx, domainBox);
            BoundaryCondition::DirichletFillGhost<double>(data_i, 0, fy, domainBox);
            BoundaryCondition::ExtrapolateCorners(data_i, domainBox);
            
            if (!(data_i.box() & xBoundBox).empty())
            {
                Box validBox = data_i.box() & domainBox;
                BoxData<double, 1> boundaryValues(validBox.edge(-Point::X()));
                boundaryValues |= Stencil<double>::CellToFace(0)(data_i);
                err[nn] = max(err[nn], boundaryValues.absMax());
            }
            if (!(data_i.box() & yBoundBox).empty())
            {
                Box validBox = data_i.box() & domainBox;
                BoxData<double, 1> boundaryValues(validBox.edge(-Point::Y()));
                boundaryValues |= Stencil<double>::CellToFace(1)(data_i);
                err[nn] = max(err[nn], boundaryValues.absMax());
            }
        }
        Reduction<double, Abs> rxn;
        rxn.reduce(&err[nn], sizeof(double));
        err[nn] = rxn.fetch();

        #if PR_VERBOSE > 0
        h5.writeMBLevel(data, "DIRICHLET_GHOST_DATA_1");
        std::cout << "Error: " << err[nn] << std::endl;
        #endif

        domainSize *= 2;
        boxSize *= 2;
    }
    for (int ii = 1; ii < numIter; ii++)
    {
        double rate = log(err[ii-1]/err[ii])/log(2.0);
        #if PR_VERBOSE > 0
        std::cout << "Convergence Rate: " << rate << std::endl;
        #endif
        EXPECT_GT(rate, 3.9);
    }
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
