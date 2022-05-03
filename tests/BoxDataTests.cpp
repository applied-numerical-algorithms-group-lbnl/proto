#include <gtest/gtest.h>
#include <cmath>
#include "Proto.H"
#include "Lambdas.H"

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

TEST(BoxData, LinearInOut) {
    Box srcBox = Box::Cube(4);                        //[(0,..,0), (3,...,3)]
    Box destBox = Box::Cube(4).shift(Point::Ones());  //[(1,...,1), (4,...,4)]

    BoxData<double,3,MEMTYPE_DEFAULT,3> Src(srcBox,7.);
    BoxData<double,2,MEMTYPE_DEFAULT,2> Dest(destBox);   //Destination data is uninitialized

    Point copyShift = Point::Ones(2); //(2,...,2)
    Box srcCopyBox = Box::Cube(3);         //[(0,...,0), (2,...,2)]

    void *buffer = new double[Src.box().size()*2*2];

    // Copy data from Src into the buffer
    Src.linearOut(buffer, srcCopyBox, {{1,2},{1,2},{0,0}});

    // ... Operate on the buffer, send it in an MPI message, etc. ...

    // Copy data from buffer into Dest
    Dest.linearIn(buffer, srcCopyBox.shift(copyShift), {{0,1},{0,1},{0,0}});

    BoxData<double,2,HOST,2> host(Dest.box());
    Dest.copyTo(host);

    for (auto it : srcCopyBox.shift(copyShift))
        EXPECT_EQ(host(it),7.);
}

TEST(BoxData, Alias) {
    Box srcBox = Box::Cube(4);
    BoxData<double,1,MEMTYPE_DEFAULT,2,3> Src(srcBox,17);
    // Alias is identical to Src and points to the same data. Changing alias will change Src.
    auto Alias = alias(Src);
    // shiftedAlias points to the same buffer as Src, but the domain is shifted by (1,...,1);
    //    (e.g. shiftedAlias[Point::Ones()] == Src[Point::Zeros] will return true.)
    auto shiftedAlias = alias(Src, Point::Ones());  //shiftedAlias points to the same data, but the associated domain
    EXPECT_TRUE(shiftedAlias.isAlias(Alias));
    EXPECT_EQ(shiftedAlias.box(),srcBox.shift(Point::Ones()));
    for (auto iter : srcBox) {
      for (int ii = 0; ii < 1; ii++)
      for (int jj = 0; jj < 2; jj++)
      for (int kk = 0; kk < 3; kk++) {
        EXPECT_EQ(Alias.data(iter,ii,jj,kk),Src.data(iter,ii,jj,kk));
        EXPECT_EQ(shiftedAlias.data((iter + Point::Ones()),ii,jj,kk),Src.data(iter,ii,jj,kk));
      }
    }
}

TEST(BoxData, Slice) {
    BoxData<int,4,MEMTYPE_DEFAULT,3,2> BD(Box::Cube(7));
    BoxData<int> DB = slice(BD,3,2,1);
    EXPECT_EQ(BD.data(3,2,1),DB.data());
    BoxData<int,5> BDslice(Box::Cube(7));
    BoxData<int,3> DBslice = slice<int,5,3>(BDslice,2);
    EXPECT_EQ(BDslice.data(2),DBslice.data());
}

int srcSize = 8, dstSize = 8, value = -1;
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
