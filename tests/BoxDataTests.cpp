#include <gtest/gtest.h>
#include <cmath>
#include "Proto.H"

using namespace Proto;
using namespace std;

TEST(BoxData, DefaultConstructor) {
    BoxData<double,2,MEMTYPE_DEFAULT,3> BD;
    EXPECT_TRUE(BD.box()==Box(Point::Zeros()));
    EXPECT_TRUE(BD.size()==0);
}

TEST(BoxData, BoxConstructor) {
    Box B = Box(Point(1,2,3,4,5,6,7));
    BoxData<int,3,MEMTYPE_DEFAULT,4,5> BD(B);
    EXPECT_TRUE(BD.box()==B);
    EXPECT_EQ(BD.size(),B.size()*3*4*5);
    EXPECT_EQ(BD.linearSize(),BD.size()*sizeof(int));
}

TEST(BoxData, Initializer) {
    Box B = Box(Point(1,2,3,4,5,6,7));
    BoxData<int,3,MEMTYPE_DEFAULT,4,5> BD(B,1337);
    EXPECT_EQ(BD.max(),1337);
    EXPECT_EQ(BD.min(),1337);
}

TEST(BoxData, Reductions) {
    size_t edge = 32;
    if ((DIM*edge) % 2) edge--;
    Box B = Box::Cube(edge);
    BoxData<int> BD(B);
    std::vector<int> buf(BD.size());
    size_t idx = B.size();
    for (auto &iter : buf) {
        iter = idx--;
        if (idx%2)
            iter *= -1;
    }
    proto_memcpy<HOST,MEMTYPE_DEFAULT>(buf.data(),BD.data(),sizeof(int)*buf.size());
    EXPECT_EQ(BD.min(),buf.front());
    EXPECT_EQ(BD.max(),buf.at(1));
    EXPECT_EQ(BD.sum(),-(buf.size()/2));
    EXPECT_EQ(BD.absMax(),-buf.front());
    BD.clampVal(0,buf.size()/2);
    EXPECT_EQ(BD.min(),0);
    EXPECT_EQ(BD.max(),buf.size()/2);
    BD.setToZero();
    EXPECT_EQ(BD.sum(),0);
    BD.setVal(-5);
    EXPECT_EQ(BD.absMax(),5);
    Reduction<int,SumAbs> rxn(true);
    BD.reduce(rxn);
    EXPECT_EQ(rxn.fetch(),5*BD.size());
}

template<typename T, unsigned int C, MemType MEM = MEMTYPE_DEFAULT>
inline bool compare(const BoxData<T,C,MEM> &bd, T val) {
    bool equal = true;
    std::vector<T> buf(bd.size());
    proto_memcpy<MEM,HOST>(bd.data(),buf.data(),sizeof(T)*buf.size());
    for (auto iter : buf) 
        if (iter != val) 
            equal = false; 
    return equal;
}

TEST(BoxData, Arithmetic) {
    Box B = Box::Cube(5);
    BoxData<int> base(B,3), sum(B,2), diff(B,1), prod(B,5), div(B,4);
    base += sum;
    EXPECT_TRUE(compare<int>(base,5));
    base -= diff;
    EXPECT_TRUE(compare<int>(base,4));
    base *= prod;
    EXPECT_TRUE(compare<int>(base,20));
    base /= div;
    EXPECT_TRUE(compare<int>(base,5));
    base.setVal(3);
    EXPECT_TRUE(compare<int>(base,3));
    base += 2;
    EXPECT_TRUE(compare<int>(base,5));
    base -= 1;
    EXPECT_TRUE(compare<int>(base,4));
    base *= 5;
    EXPECT_TRUE(compare<int>(base,20));
    base /= 4;
    EXPECT_TRUE(compare<int>(base,5));
    EXPECT_EQ(base.integrate(2),std::pow(2,DIM)*base.sum());
    std::array<int,DIM> dx;
    int factor = 1;
    for (int i=0; i<DIM;) {
        dx[i] = i+1;
        factor *= ++i;
    }
    EXPECT_EQ(base.integrate(dx),factor*base.sum());
}

TEST(BoxData, Shift) {
    Box box = Box::Cube(8);
    BoxData<int> BD(box), DB(box);
    DB.shift(Point::Ones(3));
    bool comp = true;
    auto er = DB.box().begin();
    for (auto it : BD.box()) { 
        if (it + Point::Ones(3) != *er) 
            comp = false;
        er++;
    }
    EXPECT_TRUE(comp);
}

TEST(BoxData, Alias) {
    Box box = Box::Cube(6);
    BoxData<int> BD(box);
    BoxData<int> DB = alias(BD,Point::Ones());
    EXPECT_TRUE(DB.isAlias(BD));
}

TEST(BoxData, Slice) {
    BoxData<int,4,MEMTYPE_DEFAULT,3,2> BD(Box::Cube(7));
    BoxData<int> DB = slice(BD,3,2,1);
    EXPECT_TRUE(BD.data(3,2,1)==DB.data());
    BoxData<int,5> BDslice(Box::Cube(7));
    BoxData<int,3> DBslice = slice<int,5,3>(BDslice,2);
    EXPECT_TRUE(BDslice.data(2)==DBslice.data());
}

bool testCopy(BoxData<int,2,HOST> &srcData, BoxData<int,2,HOST> &dstData,
        Point dstShift) {
    Box intersect = srcData.box().shift(dstShift) & dstData.box();
    for (auto biter = dstData.box().begin(); biter.ok(); ++biter) {
        for (int c=0; c<2; c++) {
            if (intersect.contains(*biter)) {
                if (!(srcData(*biter - dstShift, c) == dstData(*biter, c)))
                    return false;
            } else {
                if (!(dstData(*biter, c) == -1))
                    return false;
            }
        }
    }
    return true;
}

TEST(BoxData, CopyTo) {
    Box left = Box::Cube(8);
    Box right = Box::Cube(8).shift(Point::Zeros());
    Box dst = left.shift(Point::Zeros()) & right;
    Box src = dst.shift(-Point::Zeros());
    BoxData<int,2,HOST> srcData(left);
    BoxData<int,2,HOST> dstData(right);
    dstData.setVal(-1);
    srcData.copyTo(dstData, left, {0,1}, Point::Zeros(), {0,1});
    EXPECT_TRUE(testCopy(srcData, dstData, Point::Zeros()));
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
