#include <gtest/gtest.h>
#include <cmath>
#include "Proto.H"
#include "Lambdas.H"

using namespace Proto;
using namespace std;

template<typename T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
bool compareBoxData(
        BoxData<T, C, MEM, D, E>& a_src,
        BoxData<T, C, MEM, D, E>& a_dst,
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
void initBoxData(BoxData<T, C, MEM, D, E>& a_data)
{
    for (auto pt : a_data.box())
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

            a_data(pt, cc, dd, ee) = val;
        }
    }
}

TEST(BoxData, DefaultConstructor) {
    BoxData<double,2,HOST,3> BD;
    EXPECT_TRUE(BD.box()==Box(Point::Zeros()));
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
            EXPECT_EQ(hostData(p, cc, dd, ee), value);
    }
#endif
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
    constexpr unsigned int C = 2;
    int domainSize = 64;
    double dx = 1.0/domainSize;
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
    
    for (auto pt : domainBox)
    {
        for (int cc = 0; cc < C; cc++)
        {
            if (copyBox.contains(pt))
            {
                double diff = hostDst(pt + copyShift, cc) - hostSrc(pt, cc);
                EXPECT_TRUE(diff < 1e-12);
            }
        }
    }
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


TEST(BoxData, CopyToHostToHost)
{
    HDF5Handler h5;
    constexpr unsigned int COMPS = 2;
    int domainSize = 64;
    double dx = 1.0/domainSize;
    double offset = 0.125;
    Point shift = Point::Ones(7);
    double initValue = 7;
    Box srcBox = Box::Cube(domainSize);
    Box dstBoxS = srcBox.shift(shift).grow(-2);
    Box dstBoxL = srcBox.shift(shift).grow(+2);
    Box cpySrcBox = srcBox.grow(-1);
    BoxData<double, COMPS, HOST> hostSrc(srcBox);
    BoxData<double, COMPS, HOST> hostDstS(dstBoxS);
    BoxData<double, COMPS, HOST> hostDstL(dstBoxL);
    forallInPlace_p(f_phi, hostSrc, dx, offset);
    hostDstS.setVal(initValue);
    hostDstL.setVal(initValue);
    hostSrc.copyTo(hostDstL, cpySrcBox, shift);
    hostSrc.copyTo(hostDstS, cpySrcBox, shift);
    h5.writePatch(dx, hostSrc, "SRC");
    h5.writePatch(dx, hostDstL, "DST_L");
    h5.writePatch(dx, hostDstS, "DST_S");

    EXPECT_TRUE(compareBoxData(hostSrc, hostDstL, initValue, cpySrcBox, shift));
    EXPECT_TRUE(compareBoxData(hostSrc, hostDstS, initValue, cpySrcBox, shift));
}
#ifdef PROTO_CUDA
TEST(BoxData, CopyToDeviceToHost)
{
    HDF5Handler h5;
    constexpr unsigned int COMPS = 2;
    int domainSize = 64;
    double dx = 1.0/domainSize;
    double offset = 0.125;
    double initValue = 7;
    Box srcBox = Box::Cube(domainSize);
    BoxData<double, COMPS, HOST> hostSrc(srcBox);
    BoxData<double, COMPS, DEVICE> deviSrc(srcBox);
    BoxData<double, COMPS, HOST> hostDstS(srcBox.grow(-1));
    BoxData<double, COMPS, HOST> hostDstL(srcBox.grow(+1));
    forallInPlace_p(f_phi, hostSrc, dx, offset);
    forallInPlace_p(f_phi, deviSrc, dx, offset);
    hostDstS.setVal(initValue);
    hostDstL.setVal(initValue);
    deviSrc.copyTo(hostDstL);
    deviSrc.copyTo(hostDstS);
    EXPECT_TRUE(compareBoxData(hostSrc, hostDstL, srcBox, initValue));
    EXPECT_TRUE(compareBoxData(hostSrc, hostDstS, srcBox, initValue));
}

TEST(BoxData, CopyToHostToDevice)
{
    HDF5Handler h5;
    constexpr unsigned int COMPS = 2;
    int domainSize = 64;
    double dx = 1.0/domainSize;
    double offset = 0.125;
    double initValue = 7;
    Box srcBox = Box::Cube(domainSize);
    BoxData<double, COMPS, HOST> hostSrc(srcBox);
    BoxData<double, COMPS, HOST> hostDstS(srcBox.grow(-1));
    BoxData<double, COMPS, HOST> hostDstL(srcBox.grow(+1));
    BoxData<double, COMPS, DEVICE> deviDstS(srcBox.grow(-1));
    BoxData<double, COMPS, DEVICE> deviDstL(srcBox.grow(+1));
    forallInPlace_p(f_phi, hostSrc, dx, offset);
    deviDstS.setVal(initValue);
    deviDstL.setVal(initValue);
    hostSrc.copyTo(deviDstL);
    hostSrc.copyTo(deviDstS);
    proto_memcpy<DEVICE, HOST>(deviDstS.data(), hostDstS.data(), deviDstS.linearSize());
    proto_memcpy<DEVICE, HOST>(deviDstL.data(), hostDstL.data(), deviDstL.linearSize());
    EXPECT_TRUE(compareBoxData(hostSrc, hostDstL, srcBox, initValue));
    EXPECT_TRUE(compareBoxData(hostSrc, hostDstS, srcBox, initValue));
}

TEST(BoxData, CopyDeviceToDevice)
{
    HDF5Handler h5;
    constexpr unsigned int COMPS = 2;
    int domainSize = 64;
    double dx = 1.0/domainSize;
    double offset = 0.125;
    double initValue = 7;
    Box srcBox = Box::Cube(domainSize);
    BoxData<double, COMPS, HOST> hostSrc(srcBox);
    BoxData<double, COMPS, DEVICE> deviSrc(srcBox);
    BoxData<double, COMPS, HOST> hostDstS(srcBox.grow(-1));
    BoxData<double, COMPS, HOST> hostDstL(srcBox.grow(+1));
    BoxData<double, COMPS, DEVICE> deviDstS(srcBox.grow(-1));
    BoxData<double, COMPS, DEVICE> deviDstL(srcBox.grow(+1));
    forallInPlace_p(f_phi, hostSrc, dx, offset);
    forallInPlace_p(f_phi, deviSrc, dx, offset);
    deviDstS.setVal(initValue);
    deviDstL.setVal(initValue);
    hostDstS.setVal(initValue);
    hostDstL.setVal(initValue);
    deviSrc.copyTo(deviDstL);
    deviSrc.copyTo(deviDstS);
    proto_memcpy<DEVICE, HOST>(deviDstS.data(), hostDstS.data(), deviDstS.linearSize());
    proto_memcpy<DEVICE, HOST>(deviDstL.data(), hostDstL.data(), deviDstL.linearSize());
    EXPECT_TRUE(compareBoxData(hostSrc, hostDstL, srcBox, initValue));
    EXPECT_TRUE(compareBoxData(hostSrc, hostDstS, srcBox, initValue));
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
