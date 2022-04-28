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

#define COMPS 2
int value = -1;

PROTO_KERNEL_START
void f_rampHostF(const Point& a_pt, Var<int,COMPS,HOST>& a_data)
{
    Point x = a_pt + Point::Ones();
    for (int comp = 0; comp < COMPS; comp++)
    {
        a_data(comp) = (comp+1)*10 + x[0];
#if DIM > 1
        a_data(comp) = (comp+1)*100 + 10*x[0] + x[1];
#endif
#if DIM > 2
        a_data(comp) = (comp+1)*1000 + 100*x[0] + 10*x[1] + x[2];
#endif
    }
}
PROTO_KERNEL_END(f_rampHostF, f_rampHost);

PROTO_KERNEL_START
void
f_rampDeviF(const Point& a_pt, Var<int,COMPS,DEVICE>& a_data)
{
    Point x = a_pt + Point::Ones();
    for (int comp = 0; comp < COMPS; comp++)
    {
        a_data(comp) = (comp+1)*10 + x[0];
#if DIM > 1
        a_data(comp) = (comp+1)*100 + 10*x[0] + x[1];
#endif
#if DIM > 2
        a_data(comp) = (comp+1)*1000 + 100*x[0] + 10*x[1] + x[2];
#endif
    }
}
PROTO_KERNEL_END(f_rampDeviF, f_rampDevi);

int srcSize = 8, dstSize = 8;
Point shift = Point::Zeros();
//Box dst = left.shift(shift) & right;
//Box src = dst.shift(-shift);
CInterval comps = {0,COMPS-1};

TEST(BoxData, HostCopyToHost) {
    Box left = Box::Cube(srcSize);
    Box right = Box::Cube(dstSize).shift(shift);
    BoxData<int,COMPS,HOST>     srcData_h(left);
    BoxData<int,COMPS,HOST>     dstData_h(right);
    forallInPlace_p(f_rampHost, srcData_h);
    dstData_h.setVal(value);
    srcData_h.copyTo(dstData_h, left, comps, shift, comps);
    Box intersect = srcData_h.box().shift(shift) & dstData_h.box();
    for (int c=0; c<COMPS; c++) 
        for (auto biter : intersect) 
            EXPECT_EQ(srcData_h(biter-shift, c),dstData_h(biter, c));
}

#ifdef PROTO_CUDA
TEST(BoxData, HostCopyToDevice) {
    Box left = Box::Cube(srcSize);
    Box right = Box::Cube(dstSize).shift(shift);
    BoxData<int,COMPS,HOST>     srcData_h(left);
    BoxData<int,COMPS,HOST>     dstData_h(right);
    BoxData<int,COMPS,HOST>     dstData_h0(right);
    BoxData<int,COMPS,DEVICE>   dstData_d(right);
    forallInPlace_p(f_rampHost, srcData_h);
    dstData_d.setVal(value);
    dstData_h.setVal(value);
    dstData_h0.setVal(value);
    protoDeviceSynchronize(MEMTYPE_DEFAULT);
    srcData_h.copyTo(dstData_d, left, comps, shift, comps);
    protoDeviceSynchronize(MEMTYPE_DEFAULT);
    int bufferSize = dstData_d.linearSize();
    proto_memcpy<DEVICE,HOST>(dstData_d.data(), dstData_h.data(), bufferSize);
    protoDeviceSynchronize(MEMTYPE_DEFAULT);
    for (int cc = 0; cc < COMPS; cc++)
        for (auto biter : right)
            EXPECT_EQ(dstData_h(biter,cc),srcData_h(biter,cc));
}

TEST(BoxData, DeviceCopyToHost) {
    Box left = Box::Cube(srcSize);
    Box right = Box::Cube(dstSize).shift(shift);
    BoxData<int,COMPS,HOST>     srcData_h(left);
    BoxData<int,COMPS,HOST>     dstData_h(right);
    BoxData<int,COMPS,HOST>     dstData_h0(right);
    BoxData<int,COMPS,DEVICE>   srcData_d(left);
    forallInPlace_p(f_rampHost, srcData_h);
    forallInPlace_p(f_rampDevi, srcData_d);
    dstData_h.setVal(value);
    dstData_h0.setVal(value);
    protoDeviceSynchronize(MEMTYPE_DEFAULT);
    srcData_h.copyTo(dstData_h0, left, comps, shift, comps);
    srcData_d.copyTo(dstData_h, left, comps, shift, comps);
    protoDeviceSynchronize(MEMTYPE_DEFAULT);
    for (int cc = 0; cc < COMPS; cc++)
        for (auto biter : right)
            EXPECT_EQ(dstData_h(biter,cc),dstData_h0(biter,cc));
}

TEST(BoxData, DeviceCopyToDevice) {
    Box left = Box::Cube(srcSize);
    Box right = Box::Cube(dstSize).shift(shift);
    BoxData<int,COMPS,HOST>     srcData_h(left);
    BoxData<int,COMPS,HOST>     dstData_h(right);
    BoxData<int,COMPS,HOST>     dstData_h0(right);
    BoxData<int,COMPS,DEVICE>   srcData_d(left);
    BoxData<int,COMPS,DEVICE>   dstData_d(right);
    forallInPlace_p(f_rampHost, srcData_h);
    forallInPlace_p(f_rampDevi, srcData_d);
    dstData_d.setVal(value);
    dstData_h.setVal(value);
    dstData_h0.setVal(value);
    protoDeviceSynchronize(MEMTYPE_DEFAULT);
    srcData_h.copyTo(dstData_h0, left, comps, shift, comps);
    srcData_d.copyTo(dstData_d, left, comps, shift, comps);
    protoDeviceSynchronize(MEMTYPE_DEFAULT);
    dstData_d.copyTo(dstData_h);
    protoDeviceSynchronize(MEMTYPE_DEFAULT);
    for (int cc = 0; cc < COMPS; cc++)
        for (auto biter : right)
            EXPECT_EQ(dstData_h(biter,cc),dstData_h0(biter,cc));
}
#endif

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
