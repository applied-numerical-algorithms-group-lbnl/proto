#include <gtest/gtest.h>
#include "Proto.H"
#include "Lambdas.H"

PROTO_KERNEL_START
template <MemType MEM>
void consToPrim_temp(Var<double,DIM+2,MEM>& W, 
                     const Var<double, DIM+2,MEM>& U,
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

template <MemType MEM>
void forall_wrapper(Box box) {
  const double gamma = 1.4;  
  BoxData<double,DIM+2,HOST> U(Box::Cube(4),1);
  BoxData<double,DIM+2,HOST> W = forall<double,DIM+2,HOST>(consToPrim<HOST>,box,U,gamma);
  if (box==U.box())
      EXPECT_EQ(U.box(),W.box());
  else
      EXPECT_EQ(box,W.box());
  if (MEM==DEVICE) {
    BoxData<double,DIM+2,DEVICE> U_dev(U.box(),1);
    BoxData<double,DIM+2,DEVICE> W_dev = forall<double,DIM+2,DEVICE>(consToPrim,box,U_dev,gamma);
    U_dev.copyTo(U);
    W_dev.copyTo(W);
  }
  consToPrimCheck(U,W,gamma,U.box());
}

TEST(ForAll, forall) {
  Box srcBox = Box::Cube(4);
  Box destBox = Box::Cube(3);
  forall_wrapper<HOST>(srcBox);
  forall_wrapper<HOST>(destBox);
#ifdef PROTO_CUDA
  forall_wrapper<DEVICE>(srcBox);
  forall_wrapper<DEVICE>(destBox);
#endif
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

template <MemType MEM>
void forallp_wrapper(Box box) {
    const double phase = M_PI/4.;
    BoxData<double,DIM,HOST> X(Box::Cube(4),1);
    BoxData<double,1,HOST> Y = forall_p<double,1,HOST>(sineFunc,box,X,phase);
    if (box==X.box())
        EXPECT_EQ(X.box(),Y.box());
    else
        EXPECT_EQ(box,Y.box());
    if (MEM==DEVICE) {
        BoxData<double,DIM,DEVICE> X_dev(X.box(),1);
        BoxData<double,DIM,DEVICE> Y_dev = forall<double,DIM+2,DEVICE>(sineFunc,box,X_dev,gamma);
        X_dev.copyTo(X);
        Y_dev.copyTo(Y);
    }
    sineFuncCheck(X,Y,phase,Y.box());
}

TEST(ForAll, forall_p) {
  Box srcBox = Box::Cube(4);
  Box destBox = Box::Cube(3);
  forallp_wrapper<HOST>(srcBox);
  forallp_wrapper<HOST>(destBox);
#ifdef PROTO_CUDA
  forallp_wrapper<DEVICE>(srcBox);
  forallp_wrapper<DEVICE>(destBox);
#endif
}

template <MemType MEM>
void forallinplace_wrapper(Box box) {
    Box srcBox = Box::Cube(4);
    Box destBox = Box::Cube(3);
    BoxData<double,DIM+2,HOST> U(srcBox,1);
    BoxData<double,DIM+2,HOST> W(destBox);
    const double gamma = 1.4;
    forallInPlace(consToPrim,box,W,U,gamma);
    Box intersect = box & W.box();
    if (MEM==DEVICE) {
        BoxData<double,DIM+2,MEM> U_dev(srcBox,1), W_dev(destBox);
        forallInPlace(consToPrim,box,W_dev,U_dev,gamma);
        U_dev.copyTo(U,intersect);
        W_dev.copyTo(W,intersect);
    }
    consToPrimCheck(U,W,gamma,intersect);
}

TEST(ForAll, forallInPlace) {
    forallinplace_wrapper<HOST>(Box::Cube(3));
    forallinplace_wrapper<HOST>(Box::Cube(2));
#ifdef PROTO_CUDA
    forallinplace_wrapper<DEVICE>(Box::Cube(3));
    forallinplace_wrapper<DEVICE>(Box::Cube(2));
#endif
}

template <MemType MEM>
void forallinplacep_wrapper(Box box) {
    Box srcBox = Box::Cube(4);
    Box destBox = Box::Cube(3);
    BoxData<double,DIM,HOST> X(srcBox,1);
    BoxData<double,1,HOST> Y(destBox);
    const double phase = M_PI/4.;
    forallInPlace_p(sineFunc,box,Y,X,phase);
    Box intersect = box & Y.box();
    if (MEM==DEVICE) {
        BoxData<double,DIM> X_dev(srcBox,1);
        BoxData<double> Y_dev(destBox);
        forallInPlace_p(sineFunc,box,Y_dev,X_dev,phase);
        X_dev.copyTo(X,intersect);
        Y_dev.copyTo(Y,intersect);
    }
    sineFuncCheck(X,Y,phase,intersect);
}

TEST(ForAll, forallInPlace_p) {
    forallinplacep_wrapper<HOST>(Box::Cube(3));
    forallinplacep_wrapper<HOST>(Box::Cube(2));
#ifdef PROTO_CUDA
    forallinplacep_wrapper<DEVICE>(Box::Cube(3));
    forallinplacep_wrapper<DEVICE>(Box::Cube(2));
#endif
}

//TODO: Fix this test
/**
TEST(ForAll, Random) {
      const Box B0 = Box::Cube(5).shift(Point::Basis(0,-2));
      const Box B1 = Box::Cube(5);
      const Box B2 = B0 & B1;
      const Box b2 = Box::Cube(2);
      double dx = 0.1;
      auto X = forall_p<double,DIM>(iotaFunc,B0,dx);
      BoxData<int> C(B1,17);
      // forall with automatic Box

#ifdef PROTO_MEM_CHECK
      memcheck::FLUSH_CPY();
#endif
      BoxData<double,DIM> D0 = forall<double,DIM>(fooFunc,X,C);
#ifdef PROTO_MEM_CHECK
      EXPECT_EQ(memcheck::numcopies,0);
#endif
    
      EXPECT_EQ(D0.box(),B2);
      BoxData<double,DIM,HOST> xhost(B0), dhost(B2);
      BoxData<int,1,HOST> chost(B1);
      X.copyTo(xhost); C.copyTo(chost); D0.copyTo(dhost);
      for (auto it : B2)  
          EXPECT_EQ(dhost(it),chost(it)+xhost(it));
    
      // with supplied Box

#ifdef PROTO_MEM_CHECK
      memcheck::FLUSH_CPY();
#endif
      BoxData<double,DIM> D1 = forall<double,DIM>(fooFunc,b2,X,C);
#ifdef PROTO_MEM_CHECK
      EXPECT_EQ(memcheck::numcopies,0);
#endif
      EXPECT_EQ(D1.box(),b2);
      BoxData<double,DIM,HOST> host(b2);
      D1.copyTo(host);
      for (auto it : b2) 
          EXPECT_EQ(host(it),chost(it)+xhost(it));
}
*/

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
