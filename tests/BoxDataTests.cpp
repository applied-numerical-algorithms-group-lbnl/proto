#include <gtest/gtest.h>
#include <cmath>
#include "Proto.H"
#include "Lambdas.H"

using namespace Proto;
using namespace std;

template<typename T, unsigned int C, unsigned char D, unsigned char E>
bool compareBoxData(
        const BoxData<T, C, HOST, D, E>& a_src,
        const BoxData<T, C, HOST, D, E>& a_dst,
        T a_initValue,
        Box a_cpyBox,
        Point a_shift = Point::Zeros())
{
    for (auto pt : a_dst.box())
    {
        for (int ee = 0; ee < E; ee++)
        for (int dd = 0; dd < D; dd++)
        for (int cc = 0; cc < C; cc++)
        {
            if (a_cpyBox.contains(pt-a_shift))
            {
                T diff = a_dst(pt, cc, dd, ee) - a_src(pt - a_shift, cc, dd, ee);
                if (diff > 1e-12)
                {
                    std::cout << "Failure at " << pt << ", " << cc << ", " << dd << ", " << ee;
                    std::cout << " | diff: " << diff << std::endl;
                    return false;
                }
            } else {
                T diff = a_dst(pt, cc, dd, ee) - a_initValue;
                if (diff > 1e-12) { return false; }
            }
        }
    }
    return true;
}

template<typename T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
bool compareBoxData(const BoxData<T,C,MEM,D,E>& a_data, T a_val)
{
    std::vector<T> buf(a_data.size());
    proto_memcpy<MEM,HOST>(a_data.data(), buf.data(), sizeof(T)*buf.size());
    for (auto iter : buf)
    {
        if (iter != a_val)
        {
            return false;
        }
    }
    return true;
}

// This is used in lieu of forall (removes dependence on these tests on forall)
template<typename T, unsigned int C, MemType MEM, unsigned char D=1, unsigned char E=1>
void initBoxData(BoxData<T, C, MEM, D, E>& a_data)
{
    BoxData<T, C, HOST, D, E> hostData(a_data.box());
    for (auto pt : hostData.box())
    {
        for (int ee = 0; ee < E; ee++)
        for (int dd = 0; dd < D; dd++)
        for (int cc = 0; cc < C; cc++)
        {
            T val = 0;
            for (int ii = 0; ii < DIM; ii++)
            {
                val += pt[ii]*pow(100, ii);
            }
            val += (ee*100 + dd*10 + cc + 111)*pow(100, DIM);

            hostData(pt, cc, dd, ee) = val;
        }
    }
    proto_memcpy<HOST, MEM>(hostData.data(), a_data.data(), a_data.linearSize());
}

// This is used in lieu of forall (removes dependence on these tests on forall)
template<typename T, unsigned int C, MemType MEM, unsigned char D=1, unsigned char E=1>
BoxData<T, C, MEM, D, E> initBoxData(Box& a_box)
{
    BoxData<T, C, MEM, D, E> data(a_box);
    initBoxData(data);
    return data;
}

TEST(BoxData, DefaultConstructor) {
    BoxData<double,2,HOST,3> BD;
    EXPECT_TRUE(BD.box().empty());
    EXPECT_TRUE(BD.size()==0);
}

TEST(BoxData, BoxConstructor) {
    Box B = Box(Point(1,2,3,4,5,6,7));
    BoxData<int,3,HOST,4,5> BD(B);
    EXPECT_TRUE(BD.box()==B);
    EXPECT_EQ(BD.size(),B.size()*3*4*5);
    EXPECT_EQ(BD.linearSize(),BD.size()*sizeof(int));
}

template <MemType MEMTYPE>
void init_wrapper() {
    Box B = Box(Point(1,2,3,4,5,6,7));
    constexpr unsigned int C = 3;
    constexpr unsigned char D = 4;
    constexpr unsigned char E = 5;
    int value = 1337;
    BoxData<int,C,HOST,D,E> hostData(B,value);
    if (MEMTYPE==DEVICE) {
        BoxData<int,C,DEVICE,D,E> deviData(B,value);
        proto_memcpy<DEVICE, HOST>(deviData.data(), hostData.data(), deviData.linearSize());
    }
    for (auto p : B)
    {
        for (int cc = 0; cc < C; cc++) 
        for (int dd = 0; dd < D; dd++) 
        for (int ee = 0; ee < E; ee++) 
            EXPECT_EQ(hostData(p, cc, dd, ee), value);
    }
}

TEST(BoxData, Initializer) {
    init_wrapper<HOST>();
#ifdef PROTO_CUDA
    init_wrapper<DEVICE>();
#endif
}

#ifdef PROTO_MEM_CHECK
template <MemType MEMTYPE>
void move_wrapper() {
    double initValue = 1337.0;
    Box B = Box::Cube(4);

    BoxDataMemCheck::clear();
    BoxData<double, DIM, MEMTYPE> X = initBoxData<double, DIM, MEMTYPE>(B);
    EXPECT_EQ(BoxDataMemCheck::numCopies(),0);
    
    BoxDataMemCheck::clear();
    BoxData<double, DIM, MEMTYPE> Y(B, initValue);
    Y = initBoxData<double, DIM, MEMTYPE>(B);
    EXPECT_EQ(BoxDataMemCheck::numCopies(), 0);
    if (MEMTYPE==DEVICE) { 
        BoxData<double, DIM, HOST> temp(B, initValue);
        proto_memcpy<DEVICE, HOST>(Y.data(), temp.data(), Y.linearSize());
    }
    EXPECT_TRUE(compareBoxData(Y, X, initValue, B));
}

TEST(BoxData, MoveConstructor) {
    move_wrapper<HOST>();
#ifdef PROTO_CUDA
    move_wrapper<DEVICE>();
#endif
}
#endif

template <MemType MEMTYPE>
void arith_wrapper() {
    Box B = Box::Cube(5);
    std::array<int,DIM> dx;
    int factor = 1;
    for (int i=0; i<DIM;) {
        dx[i] = i+1;
        factor *= ++i;
    }

    BoxData<int, 1, MEMTYPE> base(B,3), sum(B,2), diff(B,1), prod(B,5), div(B,4);
    base += sum;
    EXPECT_TRUE(compareBoxData(base,5));
    base -= diff;
    EXPECT_TRUE(compareBoxData(base,4));
    base *= prod;
    EXPECT_TRUE(compareBoxData(base,20));
    base /= div;
    EXPECT_TRUE(compareBoxData(base,5));
    base.setVal(3);
    EXPECT_TRUE(compareBoxData(base,3));
    base += 2;
    EXPECT_TRUE(compareBoxData(base,5));
    base -= 1;
    EXPECT_TRUE(compareBoxData(base,4));
    base *= 5;
    EXPECT_TRUE(compareBoxData(base,20));
    base /= 4;
    EXPECT_TRUE(compareBoxData(base,5));
    EXPECT_EQ(base.integrate(2),std::pow(2,DIM)*base.sum());
    EXPECT_EQ(base.integrate(dx),factor*base.sum());
}

TEST(BoxData, Arithmetic) {
    arith_wrapper<HOST>();
#ifdef PROTO_CUDA
    arith_wrapper<DEVICE>();
#endif
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

TEST(BoxData, CInterval) {
    CInterval I0(1,2,3,4,5,6);
    CInterval I1{{1,2},{3,4},{5,6}};
    CInterval I2{{},{3,4},{}};
    CInterval I3{1,2};
    std::vector<CInterval> intvec{I0,I1};
    for (auto it : intvec) 
        for (int i=0; i<3; i++) {
            EXPECT_EQ(it.low(i),2*i+1);
            EXPECT_EQ(it.high(i),2*(i+1));
        } 
    EXPECT_EQ(I2.low(0),0);
    EXPECT_EQ(I2.high(0),0);
    EXPECT_EQ(I2.low(1),3);
    EXPECT_EQ(I2.high(1),4);
    EXPECT_EQ(I2.low(2),0);
    EXPECT_EQ(I2.high(2),0);
    EXPECT_EQ(I3.low(0),1);
    EXPECT_EQ(I3.high(0),2);
    EXPECT_EQ(I3.low(1),0);
    EXPECT_EQ(I3.high(1),0);
    EXPECT_EQ(I3.low(2),0);
    EXPECT_EQ(I3.high(2),0);
}

template <MemType MEMTYPE>
void linear_wrapper() {
    constexpr unsigned int C = 2;
    int domainSize = 8;
    double initValue = 7.0;
    Box domainBox = Box::Cube(domainSize);                     
    Box copyBox = Box::Cube(domainSize / 2);
    Point copyShift = Point::Ones(domainSize / 2);
    
    BoxData<double, C, MEMTYPE> src(domainBox);
    BoxData<double, C, MEMTYPE> dest(domainBox, initValue);
    initBoxData(src);

    double* buffer = (double*)proto_malloc<MEMTYPE>(copyBox.size()*C*sizeof(double));
    src.linearOut(buffer, copyBox, CInterval(0,C-1));
    dest.linearIn(buffer, copyBox.shift(copyShift), CInterval(0,C-1));
#ifdef PROTO_CUDA
        BoxData<double, C, HOST> src_host(domainBox);
        BoxData<double, C, HOST> dest_host(domainBox, initValue);
        src.copyTo(src_host);
        dest.copyTo(dest_host);
        EXPECT_TRUE(compareBoxData(src_host, dest_host, initValue, domainBox));
#else
        EXPECT_TRUE(compareBoxData(src, dest, initValue, domainBox));
#endif
}

TEST(BoxData, LinearInOut) {
    linear_wrapper<HOST>();
#ifdef PROTO_CUDA
    linear_wrapper<DEVICE>();
#endif
}

template <MemType MEMTYPE>
void alias_wrapper() {
    constexpr int C = 1;
    constexpr char D = 2;
    constexpr char E = 3;
    Box srcBox = Box::Cube(4);
    Point shift = Point::Ones();

    BoxData<double,C,MEMTYPE,D,E> hostSrc(srcBox,17);
    auto hostAlias = alias(hostSrc);
    auto hostAliasShifted = alias(hostSrc, shift);
   
    EXPECT_TRUE(hostAlias.isAlias(hostSrc));
    EXPECT_TRUE(hostAliasShifted.isAlias(hostSrc));
    EXPECT_TRUE(hostAlias.isAlias(hostAliasShifted));
    EXPECT_EQ(hostAlias.box(), srcBox);
    EXPECT_EQ(hostAliasShifted.box(), srcBox.shift(shift));
    for (auto iter : srcBox) 
        for (int ii = 0; ii < C; ii++)
        for (int jj = 0; jj < D; jj++)
        for (int kk = 0; kk < E; kk++) {
            EXPECT_EQ(hostAlias.data(iter,ii,jj,kk), hostSrc.data(iter,ii,jj,kk));
            EXPECT_EQ(hostAliasShifted.data(iter + shift,ii,jj,kk), hostSrc.data(iter,ii,jj,kk));
        }
}

TEST(BoxData, Alias) {
    alias_wrapper<HOST>();
#ifdef PROTO_CUDA 
    alias_wrapper<DEVICE>();
#endif
}

template <MemType MEMTYPE>
void slice_wrapper() {
    constexpr int  C = 4;
    constexpr char D = 3;
    constexpr char E = 2;
    
    constexpr int  c = 3;
    constexpr char d = 2;
    constexpr char e = 1;

    constexpr int cc = 2;
    int nc = 1;

    Box srcBox = Box::Cube(4);
    
    BoxData<int,C,MEMTYPE,D,E> deviSrc(srcBox);
    auto deviSlice = slice(deviSrc,c,d,e);
    EXPECT_TRUE(deviSlice.isAlias(deviSrc));

    for (auto pi : srcBox)
        EXPECT_EQ(deviSlice.data(pi,0,0,0), deviSrc.data(pi,c,d,e));

    BoxData<int,C,MEMTYPE> deviSrcV(srcBox);
    auto deviSliceV = slice<int,C,cc>(deviSrcV, nc);
    EXPECT_TRUE(deviSliceV.isAlias(deviSrcV));

    for (int ci = 0; ci < cc; ci++)
        for (auto pi : srcBox)
            EXPECT_EQ(deviSliceV.data(pi,ci), deviSrcV.data(pi,ci+nc));
}

TEST(BoxData, Slice) {
    slice_wrapper<HOST>();
#ifdef PROTO_CUDA
    slice_wrapper<DEVICE>();
#endif
}

template <MemType IN, MemType OUT>
void copyTo_wrapper() {
    constexpr unsigned int COMPS = 2;
    int domainSize = 64;
    double dx = 1.0/domainSize;
    double initValue = 7;
    Box srcBox = Box::Cube(domainSize);
    BoxData<double, COMPS, HOST> host(srcBox);
    BoxData<double, COMPS, IN> src(srcBox);
    BoxData<double, COMPS, HOST> hostDstS(srcBox.grow(-1));
    BoxData<double, COMPS, HOST> hostDstL(srcBox.grow(+1));
    BoxData<double, COMPS, OUT> destS(srcBox.grow(-1));
    BoxData<double, COMPS, OUT> destL(srcBox.grow(+1));
    initBoxData(src);
    initBoxData(host);
    hostDstS.setVal(initValue);
    hostDstL.setVal(initValue);
    if (OUT==DEVICE) {
        destS.setVal(initValue);
        destL.setVal(initValue);
        src.copyTo(destL);
        src.copyTo(destS);
        proto_memcpy<DEVICE, HOST>(destS.data(), hostDstS.data(), destS.linearSize());
        proto_memcpy<DEVICE, HOST>(destL.data(), hostDstL.data(), destL.linearSize());
    } else {
        src.copyTo(hostDstS);
        src.copyTo(hostDstL);
    }
    EXPECT_TRUE(compareBoxData(host, hostDstL, initValue, srcBox));
    EXPECT_TRUE(compareBoxData(host, hostDstS, initValue, srcBox));
}

TEST(BoxData, CopyTo)
{
    copyTo_wrapper<HOST,HOST>();
#ifdef PROTO_CUDA
    copyTo_wrapper<HOST,DEVICE>();
    copyTo_wrapper<DEVICE,HOST>();
    copyTo_wrapper<DEVICE,DEVICE>();
#endif
}

int main(int argc, char *argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
#ifdef PR_MPI
    MPI_Init(&argc, &argv);
#endif
    if (DIM > 6)
    {
        MayDay<void>::Warning("Some tests may not function for DIM > 6");
    }
#ifndef PROTO_MEM_CHECK
        MayDay<void>::Warning("PROTO_MEM_CHECK is not set. Some tests are not being run.");
#endif

    int result = RUN_ALL_TESTS();
#ifdef PR_MPI
    MPI_Finalize();
#endif
    return result;
}
