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
    constexpr unsigned int C = 3;
    constexpr unsigned char D = 4;
    constexpr unsigned char E = 5;
    int value = 1337;
    BoxData<int,C,MEMTYPE_DEFAULT,D,E> BD(B,value);
    for (auto p : B)
    {
        for (int cc = 0; cc < C; cc++)
        for (int dd = 0; dd < D; dd++)
        for (int ee = 0; ee < E; ee++)
            EXPECT_EQ(BD(p, cc, dd, ee), value);
    }
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

    BoxData<double> Src(srcBox,7.);
    BoxData<double> Dest(destBox);   //Destination data is uninitialized

    Point copyShift = Point::Ones(2); //(2,...,2)
    Box srcCopyBox = Box::Cube(3);         //[(0,...,0), (2,...,2)]

    double buffer[Src.box().size()*2*2];

    // Copy data from Src into the buffer
    Src.linearOut(buffer, srcCopyBox, CInterval(0,0));

    // ... Operate on the buffer, send it in an MPI message, etc. ...

    // Copy data from buffer into Dest
    Dest.linearIn(buffer, srcCopyBox.shift(copyShift),CInterval(0,0));

    BoxData<double,1,HOST> host(Dest.box());
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

template<typename T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
bool compareBoxData(
        BoxData<T, C, MEM, D, E>& a_src,
        BoxData<T, C, MEM, D, E>& a_dst,
        T a_initValue)
{
    for (auto pt : a_dst.box())
    {
        for (int ee = 0; ee < E; ee++)
        for (int dd = 0; dd < D; dd++)
        for (int cc = 0; cc < C; cc++)
        {
            if (a_src.box().contains(pt))
            {
                if (a_src(pt, cc, dd, ee) != a_dst(pt, cc, dd, ee)) { return false; }
            } else {
                if (a_dst(pt, cc, dd, ee) != a_initValue) { return false; }
            }
        }
    }
    return true;
}

TEST(BoxData, CopyToHostToHost)
{
    HDF5Handler h5;
    constexpr unsigned int COMPS = 2;
    int domainSize = 64;
    double dx = 1.0/domainSize;
    double offset = 0.125;
    double initValue = 7;
    Box srcBox = Box::Cube(domainSize);
    BoxData<double, 2, HOST> hostSrc(srcBox);
    BoxData<double, 2, HOST> hostDstS(srcBox.grow(-1));
    BoxData<double, 2, HOST> hostDstL(srcBox.grow(+1));
    forallInPlace_p(f_phi, hostSrc, dx, offset);
    hostDstS.setVal(initValue);
    hostDstL.setVal(initValue);
    hostSrc.copyTo(hostDstL);
    hostSrc.copyTo(hostDstS);
    EXPECT_TRUE(compareBoxData(hostSrc, hostDstL, initValue));
    EXPECT_TRUE(compareBoxData(hostSrc, hostDstS, initValue));
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
    BoxData<double, 2, HOST> hostSrc(srcBox);
    BoxData<double, 2, DEVICE> deviSrc(srcBox);
    BoxData<double, 2, HOST> hostDstS(srcBox.grow(-1));
    BoxData<double, 2, HOST> hostDstL(srcBox.grow(+1));
    forallInPlace_p(f_phi, hostSrc, dx, offset);
    forallInPlace_p(f_phi, deviSrc, dx, offset);
    hostDstS.setVal(initValue);
    hostDstL.setVal(initValue);
    deviSrc.copyTo(hostDstL);
    deviSrc.copyTo(hostDstS);
    EXPECT_TRUE(compareBoxData(hostSrc, hostDstL, initValue));
    EXPECT_TRUE(compareBoxData(hostSrc, hostDstS, initValue));
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
    BoxData<double, 2, HOST> hostSrc(srcBox);
    BoxData<double, 2, HOST> hostDstS(srcBox.grow(-1));
    BoxData<double, 2, HOST> hostDstL(srcBox.grow(+1));
    BoxData<double, 2, DEVICE> deviDstS(srcBox.grow(-1));
    BoxData<double, 2, DEVICE> deviDstL(srcBox.grow(+1));
    forallInPlace_p(f_phi, hostSrc, dx, offset);
    deviDstS.setVal(initValue);
    deviDstL.setVal(initValue);
    hostSrc.copyTo(deviDstL);
    hostSrc.copyTo(deviDstS);
    proto_memcpy<DEVICE, HOST>(deviDstS.data(), hostDstS.data(), deviDstS.linearSize());
    proto_memcpy<DEVICE, HOST>(deviDstL.data(), hostDstL.data(), deviDstL.linearSize());
    EXPECT_TRUE(compareBoxData(hostSrc, hostDstL, initValue));
    EXPECT_TRUE(compareBoxData(hostSrc, hostDstS, initValue));
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
    BoxData<double, 2, HOST> hostSrc(srcBox);
    BoxData<double, 2, DEVICE> deviSrc(srcBox);
    BoxData<double, 2, HOST> hostDstS(srcBox.grow(-1));
    BoxData<double, 2, HOST> hostDstL(srcBox.grow(+1));
    BoxData<double, 2, DEVICE> deviDstS(srcBox.grow(-1));
    BoxData<double, 2, DEVICE> deviDstL(srcBox.grow(+1));
    forallInPlace_p(f_phi, hostSrc, dx, offset);
    forallInPlace_p(f_phi, deviSrc, dx, offset);
    deviDstS.setVal(initValue);
    deviDstL.setVal(initValue);
    deviSrc.copyTo(deviDstL);
    deviSrc.copyTo(deviDstS);
    proto_memcpy<DEVICE, HOST>(deviDstS.data(), hostDstS.data(), deviDstS.linearSize());
    proto_memcpy<DEVICE, HOST>(deviDstL.data(), hostDstL.data(), deviDstL.linearSize());
    EXPECT_TRUE(compareBoxData(hostSrc, hostDstL, initValue));
    EXPECT_TRUE(compareBoxData(hostSrc, hostDstS, initValue));
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
