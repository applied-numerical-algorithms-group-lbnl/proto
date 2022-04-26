#include <gtest/gtest.h>
#include "Proto.H"

PROTO_KERNEL_START
void consToPrim_temp(Var<double,DIM+2>& W, 
                     const Var<double, DIM+2>& U,
                     double gamma)
{
  double rho = U(0);
  double v, v2=0.;
  W(0) = rho;
  
  for (int i = 1; i <= DIM; i++)
    {
      v = U(i) / rho;
      
      W(i) = v;
      v2 += v*v;
    }
  
  W(DIM+1) = (U(DIM+1) - .5 * rho * v2) * (gamma - 1.0);
}
PROTO_KERNEL_END(consToPrim_temp, consToPrim)

TEST(ForAll, forall) {
  Box srcBox = Box::Cube(4);
  BoxData<double,DIM+2> U(srcBox);
  U.setVal(1);
  const double gamma = 1.4;  
  BoxData<double,DIM+2> W = forall<double,DIM+2>(consToPrim,U,gamma);
  EXPECT_EQ(U.box(),W.box());
  BoxData<double,DIM+2,HOST> U_host(srcBox), W_host(srcBox);
  U.copyTo(U_host);
  W.copyTo(W_host);
  for (auto it : U.box()) {
    Var<double,DIM+2,HOST> u = U_host.var(it), w = W_host.var(it);
    EXPECT_EQ(u(0),w(0));
    double v2 = 0.;
    for (int i=1; i<=DIM; i++) {
      double v = u(i)/w(0);
      EXPECT_EQ(v,w(0));
      v2 += v*v;
    }
    EXPECT_EQ((u(DIM+1)-.5*u(0)*v2) * (gamma-1.),w(DIM+1));
  }
}

TEST(ForAll, forall_p) {
  
}

TEST(ForAll, forallInPlace) {
  
}

TEST(ForAll, forallInPlace_p) {
  
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
