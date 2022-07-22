#include <gtest/gtest.h>
#include <cmath>
#include "Proto.H"
#include "Lambdas.H"

using namespace Proto;
using namespace std;

template<typename T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
bool compareBoxData(
        const BoxData<T, C, MEM, D, E>& a_src,
        const BoxData<T, C, MEM, D, E>& a_dst,
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

TEST(BoxData, Initializer) {
    Box B = Box(Point(1,2,3,4,5,6,7));
    constexpr unsigned int C = 3;
    constexpr unsigned char D = 4;
    constexpr unsigned char E = 5;
    int value = 1337;
    BoxData<int,C,HOST,D,E> hostData(B,value);
    for (auto p : B)
    {
        for (int cc = 0; cc < C; cc++) 
        for (int dd = 0; dd < D; dd++) 
        for (int ee = 0; ee < E; ee++) 
            EXPECT_EQ(hostData(p, cc, dd, ee), value);
    }
#ifdef PROTO_CUDA
    BoxData<int,C,DEVICE,D,E> deviData(B,value);
    proto_memcpy<DEVICE, HOST>(deviData.data(), hostData.data(), deviData.linearSize());
    for (auto p : B)
    {
        for (int cc = 0; cc < C; cc++) 
        for (int dd = 0; dd < D; dd++) 
        for (int ee = 0; ee < E; ee++)
        {
            EXPECT_EQ(hostData(p, cc, dd, ee), value);
        }
    }
#endif
}

#ifdef PROTO_MEM_CHECK
TEST(BoxData, MoveConstructor) {
    double initValue = 1337.0;
    
    BoxDataMemCheck::clear();
    Box B = Box::Cube(4);
    BoxData<double, DIM, HOST> X_host = initBoxData<double, DIM, HOST>(B);
    EXPECT_EQ(BoxDataMemCheck::numCopies(),0);
    
    BoxDataMemCheck::clear();
    BoxData<double, DIM, HOST> Y_host(B, initValue);
    Y_host = initBoxData<double, DIM, HOST>(B);
    EXPECT_EQ(BoxDataMemCheck::numCopies(), 0);
    EXPECT_TRUE(compareBoxData(Y_host, X_host, initValue, B));

#ifdef PROTO_CUDA
    BoxDataMemCheck::clear();
    BoxData<double, DIM, DEVICE> X_devi = initBoxData<double, DIM, DEVICE>(B);
    EXPECT_EQ(BoxDataMemCheck::numCopies(),0);
    
    BoxDataMemCheck::clear();
    BoxData<double, DIM, DEVICE> Y_devi(B, initValue);
    BoxData<double, DIM, HOST>   Y_temp(B, initValue);
    Y_devi = initBoxData<double, DIM, DEVICE>(B);
    EXPECT_EQ(BoxDataMemCheck::numCopies(), 0);
    proto_memcpy<DEVICE, HOST>(Y_devi.data(), Y_temp.data(), Y_devi.linearSize());
    EXPECT_TRUE(compareBoxData(Y_temp, X_host, initValue, B));
#endif
}
#endif


TEST(BoxData, Arithmetic) {
    Box B = Box::Cube(5);
    std::array<int,DIM> dx;
    int factor = 1;
    for (int i=0; i<DIM;) {
        dx[i] = i+1;
        factor *= ++i;
    }
    
    BoxData<int, 1, HOST> base_h(B,3), sum_h(B,2), diff_h(B,1), prod_h(B,5), div_h(B,4);
    base_h += sum_h;
    EXPECT_TRUE(compareBoxData(base_h,5));
    base_h -= diff_h;
    EXPECT_TRUE(compareBoxData(base_h,4));
    base_h *= prod_h;
    EXPECT_TRUE(compareBoxData(base_h,20));
    base_h /= div_h;
    EXPECT_TRUE(compareBoxData(base_h,5));
    base_h.setVal(3);
    EXPECT_TRUE(compareBoxData(base_h,3));
    base_h += 2;
    EXPECT_TRUE(compareBoxData(base_h,5));
    base_h -= 1;
    EXPECT_TRUE(compareBoxData(base_h,4));
    base_h *= 5;
    EXPECT_TRUE(compareBoxData(base_h,20));
    base_h /= 4;
    EXPECT_TRUE(compareBoxData(base_h,5));
    EXPECT_EQ(base_h.integrate(2),std::pow(2,DIM)*base_h.sum());
    EXPECT_EQ(base_h.integrate(dx),factor*base_h.sum());

#ifdef PROTO_CUDA
    BoxData<int, 1, DEVICE> base_d(B,3), sum_d(B,2), diff_d(B,1), prod_d(B,5), div_d(B,4);
    base_d += sum_d;
    EXPECT_TRUE(compareBoxData(base_d,5));
    base_d -= diff_d;
    EXPECT_TRUE(compareBoxData(base_d,4));
    base_d *= prod_d;
    EXPECT_TRUE(compareBoxData(base_d,20));
    base_d /= div_d;
    EXPECT_TRUE(compareBoxData(base_d,5));
    base_d.setVal(3);
    EXPECT_TRUE(compareBoxData(base_d,3));
    base_d += 2;
    EXPECT_TRUE(compareBoxData(base_d,5));
    base_d -= 1;
    EXPECT_TRUE(compareBoxData(base_d,4));
    base_d *= 5;
    EXPECT_TRUE(compareBoxData(base_d,20));
    base_d /= 4;
    EXPECT_TRUE(compareBoxData(base_d,5));
    EXPECT_EQ(base_d.integrate(2),std::pow(2,DIM)*base_d.sum());
    EXPECT_EQ(base_d.integrate(dx),factor*base_d.sum());
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

TEST(BoxData, LinearInOut) {
    constexpr unsigned int C = 2;
    int domainSize = 8;
    double initValue = 7.0;
    Box domainBox = Box::Cube(domainSize);                     
    Box copyBox = Box::Cube(domainSize / 2);
    Point copyShift = Point::Ones(domainSize / 2);
    
    BoxData<double, C, HOST> hostSrc(domainBox);
    BoxData<double, C, HOST> hostDst(domainBox, initValue);
    initBoxData(hostSrc);

    double* hostBuffer = (double*)proto_malloc<HOST>(copyBox.size()*C*sizeof(double));
    hostSrc.linearOut(hostBuffer, copyBox, CInterval(0,C-1));
    hostDst.linearIn(hostBuffer, copyBox.shift(copyShift), CInterval(0,C-1));
   
    EXPECT_TRUE(compareBoxData(hostSrc, hostDst, initValue, domainBox));

#ifdef PROTO_CUDA
    BoxData<double, C, DEVICE> deviSrc(domainBox);
    BoxData<double, C, DEVICE> deviDst(domainBox, initValue);
    initBoxData(deviSrc);

    double* deviBuffer = (double*)proto_malloc<DEVICE>(copyBox.size()*C*sizeof(double));
    deviSrc.linearOut(deviBuffer, copyBox, CInterval(0,C-1));
    deviDst.linearIn( deviBuffer, copyBox.shift(copyShift), CInterval(0,C-1));
    deviDst.copyTo(hostDst);
   
    EXPECT_TRUE(compareBoxData(hostSrc, hostDst, initValue, domainBox));
#endif
}

TEST(BoxData, Alias) {
    constexpr int C = 1;
    constexpr char D = 2;
    constexpr char E = 3;
    Box srcBox = Box::Cube(4);
    Point shift = Point::Ones();

    BoxData<double,C,HOST,D,E> hostSrc(srcBox,17);
    auto hostAlias = alias(hostSrc);
    auto hostAliasShifted = alias(hostSrc, shift);
   
    EXPECT_TRUE(hostAlias.isAlias(hostSrc));
    EXPECT_TRUE(hostAliasShifted.isAlias(hostSrc));
    EXPECT_TRUE(hostAlias.isAlias(hostAliasShifted));
    EXPECT_EQ(hostAlias.box(), srcBox);
    EXPECT_EQ(hostAliasShifted.box(), srcBox.shift(shift));
    for (auto iter : srcBox) {
        for (int ii = 0; ii < C; ii++)
        for (int jj = 0; jj < D; jj++)
        for (int kk = 0; kk < E; kk++) {
            EXPECT_EQ(hostAlias.data(iter,ii,jj,kk), hostSrc.data(iter,ii,jj,kk));
            EXPECT_EQ(hostAliasShifted.data((iter + Point::Ones()),ii,jj,kk), hostSrc.data(iter,ii,jj,kk));
        }
    }

#ifdef PROTO_CUDA 
    BoxData<double,C,DEVICE,D,E> deviSrc(srcBox,17);
    auto deviAlias = alias(deviSrc);
    auto deviAliasShifted = alias(deviSrc, shift);
   
    EXPECT_TRUE(deviAlias.isAlias(deviSrc));
    EXPECT_TRUE(deviAliasShifted.isAlias(deviSrc));
    EXPECT_TRUE(deviAlias.isAlias(deviAliasShifted));
    EXPECT_EQ(deviAlias.box(), srcBox);
    EXPECT_EQ(deviAliasShifted.box(), srcBox.shift(shift));
    for (auto iter : srcBox) {
        for (int ii = 0; ii < C; ii++)
        for (int jj = 0; jj < D; jj++)
        for (int kk = 0; kk < E; kk++) {
            EXPECT_EQ(deviAlias.data(iter,ii,jj,kk), deviSrc.data(iter,ii,jj,kk));
            EXPECT_EQ(deviAliasShifted.data((iter + Point::Ones()),ii,jj,kk), deviSrc.data(iter,ii,jj,kk));
        }
    }
#endif
}

TEST(BoxData, Slice) {
    constexpr int  C = 4;
    constexpr char D = 3;
    constexpr char E = 2;
    
    constexpr int  c = 3;
    constexpr char d = 2;
    constexpr char e = 1;

    constexpr int cc = 2;
    int nc = 1;

    Box srcBox = Box::Cube(4);
    
    BoxData<int,C,HOST,D,E> hostSrc(srcBox);
    auto hostSlice = slice(hostSrc,c,d,e);
    EXPECT_TRUE(hostSlice.isAlias(hostSrc));

    for (auto pi : srcBox)
    {
        EXPECT_EQ(hostSlice.data(pi,0,0,0), hostSrc.data(pi,c,d,e));
    }

    BoxData<int,C,HOST> hostSrcV(srcBox);
    auto hostSliceV = slice<int,C,cc>(hostSrcV, nc);
    EXPECT_TRUE(hostSliceV.isAlias(hostSrcV));

    for (int ci = 0; ci < cc; ci++)
    {
        for (auto pi : srcBox)
        {
            EXPECT_EQ(hostSliceV.data(pi,ci), hostSrcV.data(pi,ci+nc));
        }
    }
    
#ifdef PROTO_CUDA
    BoxData<int,C,DEVICE,D,E> deviSrc(srcBox);
    auto deviSlice = slice(deviSrc,c,d,e);
    EXPECT_TRUE(deviSlice.isAlias(deviSrc));

    for (auto pi : srcBox)
    {
        EXPECT_EQ(deviSlice.data(pi,0,0,0), deviSrc.data(pi,c,d,e));
    }

    BoxData<int,C,DEVICE> deviSrcV(srcBox);
    auto deviSliceV = slice<int,C,cc>(deviSrcV, nc);
    EXPECT_TRUE(deviSliceV.isAlias(deviSrcV));

    for (int ci = 0; ci < cc; ci++)
    {
        for (auto pi : srcBox)
        {
            EXPECT_EQ(deviSliceV.data(pi,ci), deviSrcV.data(pi,ci+nc));
        }
    }
#endif
}


TEST(BoxData, CopyToHostToHost)
{
    constexpr unsigned int COMPS = 2;
    int domainSize = 64;
    double dx = 1.0/domainSize;
    Point shift = Point::Ones(7);
    double initValue = 7;
    Box srcBox = Box::Cube(domainSize);
    Box dstBoxS = srcBox.shift(shift).grow(-2);
    Box dstBoxL = srcBox.shift(shift).grow(+2);
    Box cpySrcBox = srcBox.grow(-1);
    BoxData<double, COMPS, HOST> hostSrc(srcBox);
    BoxData<double, COMPS, HOST> hostDstS(dstBoxS);
    BoxData<double, COMPS, HOST> hostDstL(dstBoxL);
    initBoxData(hostSrc);
    hostDstS.setVal(initValue);
    hostDstL.setVal(initValue);
    hostSrc.copyTo(hostDstL, cpySrcBox, shift);
    hostSrc.copyTo(hostDstS, cpySrcBox, shift);

    EXPECT_TRUE(compareBoxData(hostSrc, hostDstL, initValue, cpySrcBox, shift));
    EXPECT_TRUE(compareBoxData(hostSrc, hostDstS, initValue, cpySrcBox, shift));
}
#ifdef PROTO_CUDA
TEST(BoxData, CopyToDeviceToHost)
{
    constexpr unsigned int COMPS = 2;
    int domainSize = 64;
    double dx = 1.0/domainSize;
    double initValue = 7;
    Box srcBox = Box::Cube(domainSize);
    BoxData<double, COMPS, HOST> hostSrc(srcBox);
    BoxData<double, COMPS, DEVICE> deviSrc(srcBox);
    BoxData<double, COMPS, HOST> hostDstS(srcBox.grow(-1));
    BoxData<double, COMPS, HOST> hostDstL(srcBox.grow(+1));
    initBoxData(hostSrc);
    initBoxData(deviSrc);
    hostDstS.setVal(initValue);
    hostDstL.setVal(initValue);
    deviSrc.copyTo(hostDstL);
    deviSrc.copyTo(hostDstS);
    EXPECT_TRUE(compareBoxData(hostSrc, hostDstL, initValue, srcBox));
    EXPECT_TRUE(compareBoxData(hostSrc, hostDstS, initValue, srcBox));
}

TEST(BoxData, CopyToHostToDevice)
{
    constexpr unsigned int COMPS = 2;
    int domainSize = 64;
    double dx = 1.0/domainSize;
    double initValue = 7;
    Box srcBox = Box::Cube(domainSize);
    BoxData<double, COMPS, HOST> hostSrc(srcBox);
    BoxData<double, COMPS, HOST> hostDstS(srcBox.grow(-1));
    BoxData<double, COMPS, HOST> hostDstL(srcBox.grow(+1));
    BoxData<double, COMPS, DEVICE> deviDstS(srcBox.grow(-1));
    BoxData<double, COMPS, DEVICE> deviDstL(srcBox.grow(+1));
    initBoxData(hostSrc);
    deviDstS.setVal(initValue);
    deviDstL.setVal(initValue);
    hostSrc.copyTo(deviDstL);
    hostSrc.copyTo(deviDstS);
    proto_memcpy<DEVICE, HOST>(deviDstS.data(), hostDstS.data(), deviDstS.linearSize());
    proto_memcpy<DEVICE, HOST>(deviDstL.data(), hostDstL.data(), deviDstL.linearSize());
    EXPECT_TRUE(compareBoxData(hostSrc, hostDstL, initValue, srcBox));
    EXPECT_TRUE(compareBoxData(hostSrc, hostDstS, initValue, srcBox));
}

TEST(BoxData, CopyDeviceToDevice)
{
    constexpr unsigned int COMPS = 2;
    int domainSize = 64;
    double dx = 1.0/domainSize;
    double initValue = 7;
    Box srcBox = Box::Cube(domainSize);
    BoxData<double, COMPS, HOST> hostSrc(srcBox);
    BoxData<double, COMPS, DEVICE> deviSrc(srcBox);
    BoxData<double, COMPS, HOST> hostDstS(srcBox.grow(-1));
    BoxData<double, COMPS, HOST> hostDstL(srcBox.grow(+1));
    BoxData<double, COMPS, DEVICE> deviDstS(srcBox.grow(-1));
    BoxData<double, COMPS, DEVICE> deviDstL(srcBox.grow(+1));
    initBoxData(hostSrc);
    initBoxData(deviSrc);
    deviDstS.setVal(initValue);
    deviDstL.setVal(initValue);
    hostDstS.setVal(initValue);
    hostDstL.setVal(initValue);
    deviSrc.copyTo(deviDstL);
    deviSrc.copyTo(deviDstS);
    proto_memcpy<DEVICE, HOST>(deviDstS.data(), hostDstS.data(), deviDstS.linearSize());
    proto_memcpy<DEVICE, HOST>(deviDstL.data(), hostDstL.data(), deviDstL.linearSize());
    EXPECT_TRUE(compareBoxData(hostSrc, hostDstL, initValue, srcBox));
    EXPECT_TRUE(compareBoxData(hostSrc, hostDstS, initValue, srcBox));
}
#endif

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
