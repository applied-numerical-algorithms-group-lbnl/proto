#include <gtest/gtest.h>
#include "Proto.H"

PROTO_KERNEL_START
void lambdaF(Point& pt, Var<double>& data)
{
    data(0) = 0.0;
    for (int ii = 0; ii < DIM; ii++)
    {
        data(0) += pt[ii];
    }
}
PROTO_KERNEL_END(lambdaF, lambda)

TEST(InterpStencil, PiecewiseConstant) {
    int ratio = 2;
    InterpStencil<double> piecewiseConstant(ratio);
    Box iterBox = Box::Cube(ratio);  //[(0,...,0), (1,...,1)]
    for (auto destShift : iterBox)
        // The InterpStencil is indexed into using the destShift value
        piecewiseConstant(destShift) = 1.0*Shift::Zeros();

    //This function sets all the shifts / ratios automatically and effectively makes the InterpStencil read-only
    //  Using InterpStencil::operator() on a closed InterpStencil will result in an error. Use InterpStencil::get(Point destShift) for read-only access
    piecewiseConstant.close();

    Box srcBox = Box::Cube(8);
    Box computeBox = srcBox;

    auto Src = forall_p<double>(lambda,srcBox);

    BoxData<double> Dest = piecewiseConstant(Src);
    EXPECT_EQ(Dest.size(),srcBox.size()*std::pow(ratio,DIM));
    BoxData<double,1,HOST> host(Dest.box()), comp(srcBox);
    Dest.copyTo(host); Src.copyTo(comp);
    for (auto it : host.box()) {
        Point pt = it;
        for (int i=0; i<DIM; i++)
            pt[i] /= ratio;
        EXPECT_EQ(host(it),comp(pt));
    }
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
