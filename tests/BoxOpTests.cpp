#include <gtest/gtest.h>
#include "Proto.H"
#include "Lambdas.H"
#include "BoxOp_TestFlux.H"
#include "BoxOp_TestSource.H"
#include "BoxOp_TestDiag.H"

using namespace Proto;

TEST(BoxOp, Flux) {
    HDF5Handler h5;
    
    typedef BoxOp_TestFlux<double> OP;

    int domainSize = 256;
    double dx = 1.0/domainSize;
    Point k{1,2,3,4,5,6};
    Point d{1,2,3,4,5,6};
    double fluxScale = 31.0;
    double opScale = 17.0;
    OP op(dx);
    op.setFluxScale(fluxScale);

    Box B0 = Box::Cube(domainSize);
    Box B1 = B0.grow(op.ghost());
    Box B2 = B0.grow(op.auxGhost());
    BoxData<double, OP::numState()> hostSrc(B1);
    BoxData<double, OP::numAux()>   hostAux(B2);
    BoxData<double, OP::numState()> hostDst(B0);
    BoxData<double, OP::numState()> hostSln(B0);
    BoxData<double, OP::numState()> hostErr(B0);
    Array<BoxData<double, OP::numState()>, DIM> hostFlx;
    Array<BoxData<double, OP::numState()>, DIM> hostFlxSln;
    Array<BoxData<double, OP::numState()>, DIM> hostFlxErr;
    for (int dir = 0; dir < DIM; dir++)
    {
        Box fluxBox = B0.grow((Centering)dir);
        hostFlx[dir].define(fluxBox);
        hostFlxSln[dir].define(fluxBox);
        hostFlxErr[dir].define(fluxBox);
    }
    
    forallInPlace_p(f_phi, hostSrc, dx, k, d);
    hostAux.setVal(7.0);

    hostSln.setVal(0);
    for (int dir = 0; dir < DIM; dir++)
    {
        auto& F = hostFlxSln[dir];
        auto D0 = 1.0*Shift::Zeros() - 1.0*Shift::Basis(dir, -1);
        auto D1 = 1.0*Shift::Basis(dir, 1) - 1.0*Shift::Zeros();
        F |= D0(hostSrc, 1.0/dx);
        F *= hostAux;
        hostSln += D1(F, 1.0/dx);
    }
    hostSln *= (fluxScale*opScale);

    op(hostDst, hostFlx, hostSrc, hostAux, opScale);

    hostDst.copyTo(hostErr);
    hostErr -= hostSln;
    EXPECT_LT(hostErr.absMax(), 1e-12);
    for (int dir = 0; dir < DIM; dir++)
    {
        hostFlx[dir].copyTo(hostFlxErr[dir]);
        hostFlxErr[dir] -= hostFlxSln[dir];
        EXPECT_LT(hostFlxErr[dir].absMax(), 1e-12);
    }
}
TEST(BoxOp, Source) {
    typedef BoxOp_TestSource<double> OP;

    int domainSize = 256;
    double dx = 1.0/domainSize;
    Point k{1,2,3,4,5,6};
    Point d0{1,2,3,4,5,6};
    Point d1 = -d0;
    double fluxScale = 31.0;
    double opScale = 17.0;
    OP op(dx);
    op.setFluxScale(fluxScale);

    Box B0 = Box::Cube(domainSize);
    Box B1 = B0.grow(op.ghost());
    Box B2 = B0.grow(op.auxGhost());
    BoxData<double, OP::numState()> hostSrc(B1);
    BoxData<double, OP::numAux()>   hostAux(B2);
    BoxData<double, OP::numState()> hostDst(B0);
    BoxData<double, OP::numState()> hostSln(B0);
    BoxData<double, OP::numState()> hostErr(B0);
    
    forallInPlace_p(f_phi, hostSrc, dx, k, d0);
    forallInPlace_p(f_phi, hostAux, dx, k, d1);

    hostSrc.copyTo(hostSln);
    hostSln *= hostAux;
    hostSln *= opScale;
    hostSln *= fluxScale;

    op(hostDst, hostSrc, hostAux, opScale);

    hostDst.copyTo(hostErr);
    hostErr -= hostSln;
    EXPECT_LT(hostErr.absMax(), 1e-12);
}
TEST(BoxOp, Diag) {
    typedef BoxOp_TestDiag<double> OP;

    int domainSize = 256;
    double dx = 1.0/domainSize;
    Point k{1,2,3,4,5,6};
    Point d0{1,2,3,4,5,6};
    Point d1 = -d0;
    double fluxScale = 31.0;
    double opScale = 17.0;
    OP op(dx);
    op.setDiagScale(fluxScale);

    Box B0 = Box::Cube(domainSize);
    Box B1 = B0.grow(op.ghost());
    Box B2 = B0.grow(op.auxGhost());
    BoxData<double, OP::numState()> hostSrc(B1);
    BoxData<double, OP::numAux()>   hostAux(B2);
    BoxData<double, OP::numState()> hostDst(B0);
    BoxData<double, OP::numState()> hostSln(B0);
    BoxData<double, OP::numState()> hostErr(B0);
    
    forallInPlace_p(f_phi, hostSrc, dx, k, d0);
    forallInPlace_p(f_phi, hostAux, dx, k, d1);

    hostSrc.copyTo(hostSln);
    hostSln *= hostAux;
    hostSln *= opScale;
    hostSln *= fluxScale;

    op(hostDst, hostSrc, hostAux, opScale);

    hostDst.copyTo(hostErr);
    hostErr -= hostSln;
    EXPECT_LT(hostErr.absMax(), 1e-12);
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
