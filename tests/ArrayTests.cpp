#include <gtest/gtest.h>
#include "Proto.H"

using namespace Proto;

TEST(ArrayTests, Norm) {
    Array<int,2> norm = {3,4};
    EXPECT_EQ(norm.norm(),5);
}

TEST(ArrayTests, Constructors) {
    Array<double, 2> right(0);
    std::array<double,2> arr = {9.,16.};
    Array<double,2> third(arr);
    EXPECT_EQ(third,arr);
    Array<double,2> fourth(right);
    EXPECT_EQ(fourth,right);
}

TEST(ArrayTests, Modifiers) {
    Array<int,3> left, right(5);
    left.fill(5);
    EXPECT_EQ(left,right);
    right.reset();
    EXPECT_EQ(right.sum(),0);
}

TEST(ArrayTests, Arithmetic) {
    Array<int,4> left{1,2,3,4}, right = {5,6,7,8};
    Array<int,4> sum{6,8,10,12}, diff(4), prod{5,12,21,32}, div{5,3,2,2};
    EXPECT_EQ(left+right,sum);
    EXPECT_EQ(right-left,diff);
    EXPECT_EQ(left*right,prod);
    EXPECT_EQ(right/left,div);
    EXPECT_EQ(-left,-1*left);
    Array<int,4> copy(left), dummy(-4);
    left+=right;
    EXPECT_EQ(left,sum);
    left=copy;
    left-=right;
    EXPECT_EQ(left,dummy);
    left=copy;
    left*=right;
    EXPECT_EQ(left,prod);
    left=copy;
    dummy.fill(0);
    left/=right;
    EXPECT_EQ(left,dummy);
    dummy.fill(4);
    EXPECT_EQ(dummy-2,dummy/2);
    EXPECT_EQ(dummy+4,dummy*2);
    copy=dummy; 
    dummy-=2; copy/=2;
    EXPECT_EQ(copy,dummy);
    dummy+=2; copy*=2;
    EXPECT_EQ(copy,dummy);
}

TEST(ArrayTests, Modification)
{
    Array<int,6> A{1,2,3,4,5,6};
    for (int ii = 0; ii < DIM; ii++)
    {
        A[ii] = 7; 
    }
    for (int ii = 0; ii < DIM; ii++)
    {
        EXPECT_EQ(A[ii], 7);
    }
}

/* Iterator has been removed due to incorrect implementation
TEST(ArrayTests, Iterator)
{
    Array<int,6> A{1,2,3,4,5,6};
    for (auto elt : A)
    {
        elt = 7;
    }
    for (int elt : A)
    {
        EXPECT_EQ(elt, 7);
    }
}
*/
TEST(ArrayTests, Reductions) {
    Array<double,6> pi{3,-1,4,-1,5,-9}, e{-2,7,-1,8,-2,8};
    EXPECT_EQ(pi.sum(),1);
    EXPECT_EQ(e.max(),8);
    EXPECT_EQ(e.min(),-2);
    EXPECT_NE(pi.absMax(),pi.max());
    EXPECT_EQ(pi.norm(),sqrt(133));
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
