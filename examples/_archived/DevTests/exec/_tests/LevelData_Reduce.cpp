#include "Proto.H"
#include "InputParser.H"
#include "TestFunc.H"

using namespace Proto;

int main(int argc, char** argv)
{
    #ifdef PR_MPI
    MPI_Init(&argc, &argv);
    #endif
    std::cout << "DIM = " << DIM << std::endl;
    HDF5Handler h5;

    int domainSize = 64;
    double physDomainSize = 1.0;
    double rampOffset = physDomainSize / 2.0;
    int boxSize = 16;
    int ghostSize = 1;
    int numIter = 3;
    int testNum = 0;
    std::array<bool, DIM> periodicity;
    for (int dir = 0; dir < DIM; dir++) { periodicity[dir] = true; }

    InputArgs args;
    args.add("domainSize",      domainSize);
    args.add("physDomainSize",  physDomainSize);
    args.add("rampOffset",      rampOffset);
    args.add("boxSize",         boxSize);
    args.add("ghostSize",         ghostSize);
    args.add("numIter",         numIter);
    args.add("periodic_x",      periodicity[0]);
    args.add("periodic_y",      periodicity[1]);
#if DIM > 2
    args.add("periodic_z",      periodicity[2]);
#endif
    args.add("testNum",         testNum);
    args.parse(argc, argv);
    args.print();

    Point boxSizeV = Point::Ones(boxSize);
    Box domainBox = Box::Cube(domainSize);
    ProblemDomain domain(domainBox, periodicity);
    DisjointBoxLayout layout(domain, boxSizeV);

    double dx = physDomainSize / domainSize;
    if (testNum == 0)
    {
        if (procID() == 0) 
        {
            std::cout << "TEST 0: MEMTYPE = HOST" << std::endl;
        }
        LevelBoxData<double, NUMCOMPS, HOST> data0(layout, Point::Zeros());
        LevelBoxData<double, NUMCOMPS, HOST> data1(layout, Point::Ones(ghostSize));
        data0.initialize(sinProd_avg, dx); 
        data1.initialize(sinProd_avg, dx); 

        h5.writeLevel(dx, data0, "DATA_0");
        h5.writeLevel(dx, data1, "DATA_1");
        
        double absMax_0 = data0.reduce<Abs>();
        double max_0 = data0.reduce<Max>();
        double min_0 = data0.reduce<Min>();
        double sum_0 = data0.reduce<Sum>();
        
        double absMax_1 = data1.reduce<Abs>();
        double max_1 = data1.reduce<Max>();
        double min_1 = data1.reduce<Min>();
        double sum_1 = data1.reduce<Sum>();

        if (procID() == 0)
        {
            std::cout << "AbsMax (no ghost): " << absMax_0 << std::endl;
            std::cout << "Max (no ghost): " << max_0 << std::endl;
            std::cout << "Min (no ghost): " << min_0 << std::endl;
            std::cout << "Sum (no ghost): " << sum_0 << std::endl;
            
            std::cout << "AbsMax (ghost): " << absMax_1 << std::endl;
            std::cout << "Max (ghost): " << max_1 << std::endl;
            std::cout << "Min (ghost): " << min_1 << std::endl;
            std::cout << "Sum (ghost): " << sum_1 << std::endl;
        }
    } else if (testNum == 1)
    {
        if (procID() == 1) 
        {
            std::cout << "TEST 1: MEMTYPE = DEVICE" << std::endl;
        }
        LevelBoxData<double, NUMCOMPS, DEVICE> data0(layout, Point::Zeros());
        LevelBoxData<double, NUMCOMPS, DEVICE> data1(layout, Point::Ones(ghostSize));
        //data0.initialize(ramp, dx, rampOffset); 
        //data1.initialize(ramp, dx, rampOffset); 
        data0.initialize(sinProd_avg, dx); 
        data1.initialize(sinProd_avg, dx); 
        //data0.setVal(1);
        //data1.setVal(1);
        
        h5.writeLevel(dx, data0, "DATA_0");
        h5.writeLevel(dx, data1, "DATA_1");

        double absMax_0 = data0.reduce<Abs>();
        double max_0 = data0.reduce<Max>();
        double min_0 = data0.reduce<Min>();
        double sum_0 = data0.reduce<Sum>();
        
        //double absMax_1 = data1.reduce<Abs>();
        //double max_1 = data1.reduce<Max>();
        //double min_1 = data1.reduce<Min>();
        //double sum_1 = data1.reduce<Sum>();

        if (procID() == 0)
        {
            std::cout << "AbsMax (no ghost): " << absMax_0 << std::endl;
            std::cout << "Max (no ghost): " << max_0 << std::endl;
            std::cout << "Min (no ghost): " << min_0 << std::endl;
            std::cout << "Sum (no ghost): " << sum_0 << std::endl;
            
            //std::cout << "AbsMax (ghost): " << absMax_1 << std::endl;
            //std::cout << "Max (ghost): " << max_1 << std::endl;
            //std::cout << "Min (ghost): " << min_1 << std::endl;
            //std::cout << "Sum (ghost): " << sum_1 << std::endl;
        }
    }
}
