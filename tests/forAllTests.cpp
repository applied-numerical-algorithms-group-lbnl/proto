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
  
  W(DIM+1) = (U(DIM+1)-.5*rho*v2) * (gamma-1.);
}
PROTO_KERNEL_END(consToPrim_temp, consToPrim)

void consToPrimCheck(const BoxData<double,DIM+2,HOST> &U, const BoxData<double,DIM+2,HOST> &W, const double gamma, Box domain) {
  for (auto it : domain) {
    Var<double,DIM+2,HOST> u = U.var(it), w = W.var(it);
    EXPECT_EQ(u(0),w(0));
    double v2 = 0.;
    for (int i=1; i<=DIM; i++) {
      double v = u(i)/u(0);
      EXPECT_EQ(v,w(0));
      v2 += v*v;
    }
    EXPECT_EQ((u(DIM+1)-.5*u(0)*v2) * (gamma-1.),w(DIM+1));
  }
}

TEST(ForAll, forall) {
  Box srcBox = Box::Cube(4);
  BoxData<double,DIM+2> U(srcBox,1);
  const double gamma = 1.4;  
  BoxData<double,DIM+2> W = forall<double,DIM+2>(consToPrim,U,gamma);
  EXPECT_EQ(U.box(),W.box());
  BoxData<double,DIM+2,HOST> U_host(srcBox), W_host(srcBox);
  U.copyTo(U_host);
  W.copyTo(W_host);
  consToPrimCheck(U_host,W_host,gamma,srcBox);

  Box destBox = Box::Cube(3);
  BoxData<double,DIM+2> W2 = forall<double,DIM+2>(consToPrim,destBox,U,gamma);
  EXPECT_EQ(W2.box(),destBox);
  BoxData<double,DIM+2,HOST> W2_host(destBox);
  W2.copyTo(W2_host);
  consToPrimCheck(U_host,W2_host,gamma,destBox);
}

PROTO_KERNEL_START
void sineFunc_temp(Point& a_p,
                   Var<double>&     a_Y,
                   const Var<double,DIM>& a_X,
                   double           phase)
{
  a_Y(0) = 0.0;
  for (int ii = 0; ii < DIM; ii++)
  {
    a_Y(0) += sin(a_X(ii) + phase);
  }
}
PROTO_KERNEL_END(sineFunc_temp, sineFunc)

void sineFuncCheck(const BoxData<double,DIM,HOST> &X, const BoxData<double,1,HOST> &Y, double phase, const Box & domain) {
    for (auto it : domain) {
        double comp = 0.;
        Var<double,DIM,HOST> x = X.var(it);
        Var<double,1,HOST> y = Y.var(it);
        for (int i=0; i<DIM; i++) 
            comp += sin(x(i)+phase);
        EXPECT_EQ(comp,y(0));
    }
}

TEST(ForAll, forall_p) {
    Box srcBox = Box::Cube(4);
    BoxData<double,DIM> X(srcBox,1);
    const double phase = M_PI/4.;
    BoxData<double> Y = forall_p<double>(sineFunc,X,phase);
    EXPECT_EQ(X.box(),Y.box());
    BoxData<double,DIM,HOST> X_host(X.box());
    BoxData<double,1,HOST> Y_host(Y.box());
    X.copyTo(X_host);
    Y.copyTo(Y_host);
    sineFuncCheck(X_host,Y_host,phase,Y.box());
    
    Box destBox = Box::Cube(3);
    BoxData<double> Y2 = forall_p<double>(sineFunc,destBox,X,phase);
    EXPECT_EQ(destBox,Y2.box());
    BoxData<double,1,HOST> Y2_host(destBox);
    Y2.copyTo(Y2_host);
    sineFuncCheck(X_host,Y2_host,phase,Y2.box());
}

TEST(ForAll, forallInPlace) {
    Box srcBox = Box::Cube(4);
    Box destBox = Box::Cube(3);
    BoxData<double,DIM+2> U(srcBox,1);
    BoxData<double,DIM+2> W(destBox);
    const double gamma = 1.4;
    forallInPlace(consToPrim,W,U,gamma);
    BoxData<double,DIM+2,HOST> U_host(U.box()), W_host(W.box());
    U.copyTo(U_host);
    W.copyTo(W_host);
    consToPrimCheck(U_host,W_host,gamma,destBox);

    Box computeBox = Box::Cube(2);
    forallInPlace(consToPrim,computeBox,W,U,gamma);
    Box intersect = computeBox & W.box();
    W.copyTo(W_host,intersect);
    consToPrimCheck(U_host,W_host,gamma,intersect);
}

TEST(ForAll, forallInPlace_p) {
    Box srcBox = Box::Cube(4);
    Box destBox = Box::Cube(3);
    BoxData<double,DIM> X(srcBox,1);
    BoxData<double> Y(destBox);
    const double phase = M_PI/4.;
    forallInPlace_p(sineFunc,Y,X,phase);
    BoxData<double,DIM,HOST> X_host(X.box());
    BoxData<double,1,HOST> Y_host(Y.box());
    X.copyTo(X_host);
    Y.copyTo(Y_host);
    sineFuncCheck(X_host,Y_host,phase,Y.box());

    Box computeBox = Box::Cube(2);
    forallInPlace_p(sineFunc,computeBox,Y,X,phase);
    Box intersect = computeBox & Y.box();
    Y.copyTo(Y_host,intersect);
    sineFuncCheck(X_host,Y_host,phase,intersect);
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
