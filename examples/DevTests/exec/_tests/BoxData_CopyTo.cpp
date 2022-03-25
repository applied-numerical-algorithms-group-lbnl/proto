#include "Proto.H"
#include "InputParser.H"

#define NUMCOMPS 2

using namespace Proto;

PROTO_KERNEL_START
void
f_rampF(const Point& a_pt, Var<double,NUMCOMPS>& a_data)
{
    Point x = a_pt + Point::Ones();
    for (int comp = 0; comp < NUMCOMPS; comp++)
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
PROTO_KERNEL_END(f_rampF, f_ramp);

template <MemType MEM>
void testCopy(
        BoxData<double, NUMCOMPS, MEM>& a_srcData,
        BoxData<double, NUMCOMPS, MEM>& a_dstData,
        Point a_dstShift)
{
    bool pass = true;
    Box intersect = a_srcData.box().shift(a_dstShift) & a_dstData.box();
    for (auto biter = a_dstData.box().begin(); biter.ok(); ++biter)
    {
        for (int c = 0; c < NUMCOMPS; c++)
        {
            if (intersect.contains(*biter))
            {
                pass &= (a_srcData(*biter - a_dstShift, c) == a_dstData(*biter, c));
                if (!pass)
                {
                    std::cout << "TEST FAILED | dstBox point: " << *biter;
                    std::cout << " | srcData value: " << a_srcData(*biter - a_dstShift, c);
                    std::cout << " | dstData value: " << a_dstData(*biter, c);
                    std::cout << std::endl;
                    std::abort();
                }
            } else {
                pass &= (a_dstData(*biter, c) == -1);
                if (!pass)
                {
                    std::cout << "TEST FAILED | dstBox point: " << *biter;
                    std::cout << " | dstData value: " << a_dstData(*biter, c);
                    std::cout << " (should be -1)";
                    std::cout << std::endl;
                    std::abort();
                }
            }
        }
    }
    if (pass)
    {
        std::cout << "TEST PASSED" << std::endl;
    }
}

int main(int argc, char** argv)
{
    
    #ifdef PR_MPI
    MPI_Init(&argc, &argv);
    #endif

    // SETUP
    HDF5Handler h5;

    int srcBoxSize = 8;
    int dstBoxSize = 8;
    Point dstShift = Point::Zeros();
    Point cpyShift = Point::Zeros();
    int testNum = 0;

    InputArgs args;
    args.add("srcBoxSize",  srcBoxSize);
    args.add("dstBoxSize",  srcBoxSize);
    args.add("dstShift_x",  dstShift[0]);
    args.add("cpyShift_x",  cpyShift[0]);
#if DIM > 1
    args.add("dstShift_y",  dstShift[1]);
    args.add("cpyShift_y",  cpyShift[1]);
#endif
#if DIM > 2
    args.add("dstShift_z",  dstShift[2]);
    args.add("cpyShift_z",  cpyShift[2]);
#endif
    args.add("testNum", testNum);
    args.parse(argc, argv);
    args.print();

    std::array<double, DIM> dx;
    dx.fill(1);
    Box srcBox0 = Box::Cube(srcBoxSize);
    Box dstBox0 = Box::Cube(dstBoxSize).shift(dstShift);
    
    CInterval srcComps = {0, NUMCOMPS-1};
    CInterval dstComps = {0, NUMCOMPS-1};

    Box dstBox = srcBox0.shift(cpyShift) & dstBox0;
    Box srcBox = dstBox.shift(-cpyShift);

    std::cout << "Copying from " << srcBox0 << " to " << dstBox0 << " with shift " << cpyShift << std::endl;
    std::cout << "Src intersection box: " << srcBox << std::endl;
    std::cout << "Dst intersection box: " << dstBox << std::endl;
    std::cout << "Src components: " << srcComps << std::endl;
    std::cout << "Dst components: " << dstComps << std::endl;

    if (testNum == 0)
    {
        std::cout << " Test 0: Host - Host Copy" << std::endl;
        BoxData<double, NUMCOMPS, HOST> srcData(srcBox0);
        BoxData<double, NUMCOMPS, HOST> dstData(dstBox0);
        dstData.setVal(-1);
        forallInPlace_p(f_ramp, srcData);
        srcData.copyTo(dstData, srcBox0, srcComps, cpyShift, dstComps);
        testCopy(srcData, dstData, cpyShift);
    }
#ifdef PROTO_CUDA
    else if (testNum == 1)
    {
        std::cout << " Test 1: Host - Device Copy" << std::endl;
        BoxData<double, NUMCOMPS, HOST>     srcData_h(srcBox0);
        BoxData<double, NUMCOMPS, DEVICE>   dstData_d(dstBox0);
        BoxData<double, NUMCOMPS, DEVICE>   srcData_d(srcBox0);
        dstData_d.setVal(-1);
        forallInPlace_p(f_ramp, srcData_h);
        forallInPlace_p(f_ramp, srnData_d);
        srcData_h.copyTo(dstData_d, srcBox0, srcComps, cpyShift, dstComps);
        testCopy(srcData_d, dstData_d, cpyShift);
    } else if (testNum == 2)
    {
        std::cout << " Test 2: Device - Host Copy" << std::endl;
    } else if (testNum == 3)
    {
        std::cout << " Test 3: Device - Device Copy" << std::endl;
    }
#endif
    else {
        std::cout << " Unrecognized Test Number: " << testNum << std::endl;
    }
    #ifdef PR_MPI
    MPI_Finalize();
    #endif
}
