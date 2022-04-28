#include <gtest/gtest.h>
#include "Proto.H"
#include <numeric>

TEST(Stencil, Build) {
    Stencil<double> L = (-2.0*DIM)*Shift::Zeros();
    for (int dir = 0; dir < DIM; dir++)
      L += 1.0*Shift::Basis(dir,1) + 1.0*Shift::Basis(dir,-1);
    EXPECT_EQ(L,Stencil<double>::Laplacian());
}

TEST(Stencil, Utility) {
    Stencil<int> test = 1*Shift::Zeros();
    for (int dir=0; dir<DIM; dir++)
        test += 1*Shift::Basis(dir,DIM+dir) + 1*Shift::Basis(dir,dir-DIM);
    EXPECT_EQ(test.span(),Box(Point{-DIM,1-DIM,2-DIM},Point{DIM,DIM+1,DIM+2}));
    EXPECT_EQ(test.ghost(),Point(DIM,DIM+1));
    Stencil<int> L = Stencil<int>::Laplacian();
    Stencil<int> orig(L), trans(L);
    L.invert(0);
    auto ti = orig.offsets().begin();
    for (auto it : L.offsets())
        EXPECT_EQ(-it[0],(*(ti++))[0]);
    L.transpose(0,1);
    EXPECT_EQ(L,trans);
}

TEST(Stencil, Operators) {
    Stencil<int> left  =  1*Shift::Basis(0,-1)-1*Shift::Zeros()+1*Shift::Basis(0,1);
    Stencil<int> right = -1*Shift::Basis(1,-1)+1*Shift::Zeros()-1*Shift::Basis(1,1);
    Stencil<int> sum(left+right);
    EXPECT_EQ(left-right,left+(-1*right));
    Stencil<int> prod(left*right);
    EXPECT_EQ(left.size()*right.size(),prod.size());
    std::vector<int> mult;
    double add = std::accumulate(left.coefs().begin(),left.coefs().end(),0.)
        + std::accumulate(right.coefs().begin(),right.coefs().end(),0.);
    EXPECT_EQ(add,std::accumulate(sum.coefs().begin(),sum.coefs().end(),0.));
    Stencil<int> temp(sum);
    sum -= right;
    EXPECT_EQ(sum,temp+(-1*right));
    for (auto lit : left.coefs())
        for (auto rit : right.coefs())
            mult.push_back(lit*rit);
    EXPECT_EQ(prod.coefs(),mult);
    left *= right;
    EXPECT_EQ(left,prod);
    left += right;
    EXPECT_EQ(left,prod+right);
    left -= right;
    Stencil<int> scal = right*5;
    auto ti = scal.coefs().begin();
    for (auto it : right.coefs())
        EXPECT_EQ(it*5,*(ti++));
    right *= 5;
    EXPECT_EQ(scal,right);
}

PROTO_KERNEL_START
void sineFunc_temp(Point& pt, Var<double>& src, double dx) {
  src(0) = sin(pt[0]*dx); //set Src = sin(x)
}
PROTO_KERNEL_END(sineFunc_temp, sineFunc)

TEST(Stencil, Apply) {
    auto L = Stencil<double>::Laplacian(); //2nd order 2*DIM + 1 Point Laplacian
      //  Initialize some source data
      //  For this example, we define our input based on our desired ouput: an 8^DIM Point cube.
    int rangeSize = 8;
    Box rangeBox = Box::Cube(rangeSize);  // [(0,...,0), (7,...,7)]
      //  The domainBox can be computed from the desired range, tanking ghost (halo) cells into account
    Box domainBox = L.domain(rangeBox);  // L needs 1 layer of ghost cells, so the domainBox is: [ (-1,...,-1), (8,....,8)]
    Box computeBox = rangeBox.grow(-1);
   //  Initialize the data using forallV
    double dx = 0.5*M_PI/rangeSize; //assuming a range of PI/2
    BoxData<double> Src = forall_p<double>(sineFunc, domainBox, dx); //don't forget the domain Box!

    //  Apply L to Src
      //  Method 1: Build a new BoxData from the Stencil computation
    BoxData<double> Dest_0 = L(Src);
      //  When no Box input is specified, the output will be defined on the largest possible Box (e.g. rangeBox)
    BoxData<double> Dest_1 = L(Src, computeBox);
    EXPECT_EQ(Dest_0.box(),rangeBox);
    EXPECT_EQ(Dest_1.box(),computeBox);
    BoxData<double,1,HOST> host_0(Dest_0.box()), host_1(Dest_1.box());
    Dest_1.copyTo(host_1);
    Dest_0.copyTo(host_0);
    for (auto it : host_1.box())
        EXPECT_EQ(host_0(it),host_1(it));
      //  When a VALID Box input is specified, the output will be restricted to that Box.
      //  *Usually* it makes more sense to let Stencil do its magic and range Box for you.

      //  Method 2: Apply Stencil computation output to an existing BoxData
    BoxData<double> Dest_3(rangeBox);
    BoxData<double> Dest_4(rangeBox,0);
      //  REPLACE the data in the destination with the output of L(Src)
    Dest_3 |= L(Src);
      //  ADD the output of L(Src) to the existing data in the destination
    Dest_4 += L(Src);
      //  Both the ADD and REPLACE operations can specify a compute Box as well if desired
      //  Again, specifying the Box input is not recommended unless it is necessary
    EXPECT_EQ(Dest_3.box(),Dest_4.box());
    BoxData<double,1,HOST> host_3(Dest_3.box()), host_4(Dest_4.box());
    Dest_3.copyTo(host_3);
    Dest_4.copyTo(host_4);
    for (auto it : host_3.box())
        EXPECT_EQ(host_3(it),host_4(it));
}

PROTO_KERNEL_START
void initFunc_temp(Point &pt, Var<double> &src) {
  src(0) = 0.;
  for (int ii = 0; ii < DIM; ii++)
    src(0) += pt[ii];
}
PROTO_KERNEL_END(initFunc_temp, initFunc)

TEST(Stencil, SourceRefine) {
    //  This example Stencil computes a linear average from fine data source data onto a coarse grid
    Stencil<double> Avg;
    int refRatio = 2;
    Box offsetBox = Box::Cube(refRatio);
    for (auto iter : offsetBox)
      Avg += 1.*Shift(iter);
    Avg *= (1./offsetBox.size());
    //  WARNING: Stencils with different src/dest refinement ratios cannot be added, subtracted, etc.
    //    When building a Stencil with non-trivial refinement, set the refinement ratio last.
    Avg.srcRatio() = Point::Ones(refRatio);

    int rangeSize = 8;
    Box indexBox = Box::Cube(rangeSize);          //[(0,...,0), (7,....,7)]
    Box domainBox = indexBox.refine(refRatio);   //[(0,...,0), (15,...,15)]
    auto Src = forall_p<double>(initFunc,domainBox);

    //  The index box represents the points at which to compute the Stencil.
    //    When the source is refined by r, the data used to compute the solution at point p and p + 1
    //    In the index box is shifted by r.
    
    //    In this example, r = 2. In 1D, the average at Point (0) is computed by:
    //      avg(0) = (1.0*src(0) + 1.0*src(1))/2.
    //    At the next Point, the forumula is:
    //      avg(1) = (1.0*src(2) + 1.0*src(3))/2.
    //    Manifesting a shift in the source data by r = 2.

    BoxData<double> Dest_0 = Avg(Src,indexBox);
    //  OR
    BoxData<double> Dest_1 = Avg(Src); //Stencil automatically determines the largest possible Box for Dest_1, given the data available.
    //  The "|=" and "+=" operators can be used here as well with identical semantics.
    //  For an example illustrating the usage of the destination refinement ratio, see the documentation for InterpStencil
    EXPECT_EQ(indexBox,Dest_1.box());
    BoxData<double,1,HOST> a(indexBox), b(indexBox);
    Dest_0.copyTo(a); Dest_1.copyTo(b);
    for (auto pt : a.box())
        EXPECT_EQ(a(pt),b(pt));
}

TEST(Stencil, DestRefine) {
    int refRatio = 2;
    //  This is a very simple stencil; it scales a Point by 7.
    double coef = 7.;
    Stencil<double> S = coef*Shift::Zeros();
    S.destRatio() = Point::Ones(2);
    S.destShift() = Point::Ones();

    int domainSize = 8;
    Box indexBox = Box::Cube(domainSize);           //[(0,...,0), (7,....,7)]
    Box rangeBox = indexBox.refine(refRatio);      //[(0,...,0), (15,...,15)]
    BoxData<double> Src(indexBox,1.);
    BoxData<double> Dest(rangeBox,0.);

    Dest |= S(Src);
    EXPECT_EQ(Src.box().refine(2),Dest.box());
    BoxData<double,1,HOST> host(rangeBox);
    Dest.copyTo(host);
    for (auto it : Dest.box())
        if (host(it))
            EXPECT_EQ(coef,host(it));
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
