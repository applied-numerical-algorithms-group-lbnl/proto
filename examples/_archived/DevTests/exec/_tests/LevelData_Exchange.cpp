#include "Proto.H"
#include "InputParser.H"

#define NUMCOMPS 1

using namespace Proto;

PROTO_KERNEL_START
void
f_rampHostF(const Point& a_pt, Var<int,NUMCOMPS, HOST>& a_data, Box a_domainBox)
{
    Point x = (a_domainBox % a_pt) + Point::Ones();
    
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
PROTO_KERNEL_END(f_rampHostF, f_rampHost);

PROTO_KERNEL_START
void
f_rampDeviF(const Point& a_pt, Var<int,NUMCOMPS, DEVICE>& a_data, Box a_domainBox)
{
    Point x = (a_domainBox % a_pt) + Point::Ones();
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
PROTO_KERNEL_END(f_rampDeviF, f_rampDevi);

int main(int argc, char** argv)
{
#ifdef PR_MPI
    MPI_Init(&argc, &argv);
#endif

    HDF5Handler h5;

    int domainSize = 64;
    int boxSize = 16;
    int ghostSize = 1;
    std::array<bool, DIM> periodicity;
    int testNum = 0;
    for (int dir = 0; dir < DIM; dir++) { periodicity[dir] = true; }

    InputArgs args;
    args.add("domainSize",      domainSize);
    args.add("boxSize",         boxSize);
    args.add("ghostSize",       ghostSize);
    args.add("testNum",         testNum);
    args.add("periodic_x",      periodicity[0]);
    args.add("periodic_y",      periodicity[1]);
#if DIM > 2
    args.add("periodic_z",      periodicity[2]);
#endif
    args.parse(argc, argv);
    args.print();

    std::array<double, DIM> dx;
    dx.fill(1.0);

    Point boxSizeV = Point::Ones(boxSize);
    Box domainBox = Box::Cube(domainSize);
    ProblemDomain domain(domainBox, periodicity);
    DisjointBoxLayout layout(domain, boxSizeV);

    if (testNum == 0)
    {
        if (procID() == 0)
        {
            std::cout << "RUNNING TEST 0: MEMTYPE = HOST" << std::endl;
        }
        LevelBoxData<int, NUMCOMPS, HOST> srcData(layout, Point::Zeros());
        LevelBoxData<int, NUMCOMPS, HOST> dstData(layout, Point::Ones(ghostSize));
        LevelBoxData<int, NUMCOMPS, HOST> slnData(layout, Point::Ones(ghostSize));

        dstData.setVal(-1);
        srcData.initialize(f_rampHost, domainBox); 
        slnData.initialize(f_rampHost, domainBox);
        for (auto iter = layout.begin(); iter.ok(); ++iter)
        {
            auto& src_i = srcData[*iter];
            auto& dst_i = dstData[*iter];
            src_i.copyTo(dst_i);
        }
        
        dstData.exchange();
       
        bool PASS = true;
        for (auto iter = layout.begin(); iter.ok(); ++iter)
        {
            auto& dst_i = dstData[*iter];
            auto& sln_i = slnData[*iter];
            for (int cc = 0; cc < NUMCOMPS; cc++)
            {
                for (auto biter = dst_i.box().begin(); biter.ok(); ++biter)
                {
                    bool pass = (dst_i(*biter, cc) == sln_i(*biter, cc));
                    if (!pass && procID() == 0)
                    {
                        std::cout << "TEST FAILED | patch: " << (*iter).global();
                        std::cout << " | dst(" << *biter << ", " << cc << ") = " << dst_i(*biter, cc);
                        std::cout << " != sln(" << *biter << ", " << cc << ") = " << sln_i(*biter, cc) << std::endl;
                    }
                    PASS &= pass;
                }
            }
        }
        if (procID() == 0 && PASS)
        {
            std::cout << "TEST PASSED" << std::endl;
        }
        
    }
    else if (testNum == 1)
    {
        if (procID() == 0)
        {
            std::cout << "RUNNING TEST 1: MEMTYPE = DEVICE" << std::endl;
        }
#ifdef PROTO_ACCEL
        LevelBoxData<int, NUMCOMPS, DEVICE> srcData(layout, Point::Zeros());
        LevelBoxData<int, NUMCOMPS, DEVICE> dstData(layout, Point::Ones(ghostSize));
        LevelBoxData<int, NUMCOMPS, DEVICE> slnData(layout, Point::Ones(ghostSize));

        dstData.setVal(-1);
        srcData.initialize(f_rampDevi, domainBox); 
        slnData.initialize(f_rampDevi, domainBox);
        for (auto iter = layout.begin(); iter.ok(); ++iter)
        {
            auto& src_i = srcData[*iter];
            auto& dst_i = dstData[*iter];
            src_i.copyTo(dst_i);
        }
        
        dstData.exchange();
       
        bool PASS = true;
        for (auto iter = layout.begin(); iter.ok(); ++iter)
        {
            auto& dst_d = dstData[*iter];
            auto& sln_d = slnData[*iter];
            BoxData<int, NUMCOMPS, HOST> dst_h(dst_d.box());
            BoxData<int, NUMCOMPS, HOST> sln_h(sln_d.box());
            dst_d.copyTo(dst_h);
            sln_d.copyTo(sln_h);
            for (int cc = 0; cc < NUMCOMPS; cc++)
            {
                for (auto biter = dst_h.box().begin(); biter.ok(); ++biter)
                {
                    bool pass = (dst_h(*biter, cc) == sln_h(*biter, cc));
                    if (!pass)
                    {
                        pout() << "TEST FAILED | patch: " << (*iter).global();
                        pout() << " | dst(" << *biter << ", " << cc << ") = " << dst_h(*biter, cc);
                        pout() << " != sln(" << *biter << ", " << cc << ") = " << sln_h(*biter, cc) << std::endl;
                    }
                    PASS &= pass;
                }
            }
        }
        if (PASS)
        {
                std::cout << "TEST PASSED (proc " << procID() << ")" << std::endl;
        } else {
                std::cout << "TEST FAILED (proc " << procID() << ")" << std::endl;
        }
#else
        if (procID() == 0)
        {
            std::cout << "TEST NOT RUN: No device option is active." << std::endl;
        }
#endif
    }
#ifdef PR_MPI
    MPI_Finalize();
#endif
}
