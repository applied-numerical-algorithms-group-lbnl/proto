#include <gtest/gtest.h>
#include "Proto.H"
#include "Lambdas.H"

using namespace Proto;

TEST(StencilLib, CornersToCells) {
#ifdef PR_HDF5
    HDF5Handler h5;
#endif
    auto S = Stencil<double>::CornersToCells(4);

    int boxSize = 64;
    int N = 3;
    double err[N];
    for (int nn = 0; nn < N; nn++)
    {
        Box B1 = Box::Cube(boxSize);
        Box B0 = S.domain(B1);


        Array<double, DIM> dx = 1.0/boxSize;
        Array<double, DIM> offset = 0.0;
        Array<double, DIM> k = 1.0;
        auto srcData = forall_p<double>(f_phi_corner, B0, dx, k, offset); 
        auto slnData = forall_p<double>(f_phi_avg, B1, dx, k, offset); 

        BoxData<double> dstData = S(srcData);
        BoxData<double> errData(B1);
        dstData.copyTo(errData);
        errData -= slnData;

#ifdef PR_HDF5
        h5.writePatch(dx, errData, "C2C_ERR_%i", nn);
        h5.writePatch(dx, srcData, "C2C_SRC_%i", nn);
#endif

        err[nn] = errData.absMax();
        boxSize *= 2;
#if PR_VERBOSE > 0
        std::cout << "Error: " << err[nn] << std::endl;
        srcData.printData();
        dstData.printData();
        slnData.printData();
        errData.printData();
#endif
    }

    for (int ii = 1; ii < N; ii++)
    {
        double rate = log(err[ii-1]/err[ii])/log(2.0);
#if PR_VERBOSE > 0
        std::cout << "Convergence Rate: " << rate << std::endl;
#endif
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
