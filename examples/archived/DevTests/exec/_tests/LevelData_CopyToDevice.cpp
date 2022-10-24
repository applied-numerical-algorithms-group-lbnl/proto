#include "Proto.H"
#include "TestFunc.H"
#include "InputParser.H"

using namespace Proto;

int main(int argc, char** argv)
{
    #ifdef PR_MPI
    MPI_Init(&argc, &argv);
    #endif

    #ifdef PROTO_ACCEL
    // SETUP
    HDF5Handler h5;

    int domainSize = 64;
    int boxSize = 16;
    std::array<bool, DIM> periodicity;
    for (int dir = 0; dir < DIM; dir++) { periodicity[dir] = true; }

    InputArgs args;
    args.add("domainSize",      domainSize);
    args.add("boxSize",         boxSize);
    args.add("periodic_x",      periodicity[0]);
    args.add("periodic_y",      periodicity[1]);
#if DIM > 2
    args.add("periodic_z",      periodicity[2]);
#endif
    args.parse(argc, argv);
    args.print();
   
    std::array<double, DIM> dx;
    dx.fill(1.0);

    Box domainBox = Box::Cube(domainSize);
    Point boxSizeVect = Point::Ones(boxSize);
    ProblemDomain domain(domainBox, periodicity);
    DisjointBoxLayout layout(domain, boxSizeVect);
    
    pout() << "Defining device data" << std::endl;
    LevelBoxData<int, NUMCOMPS, DEVICE> data_d(layout, Point::Ones(1));
    pout() << "Defining host data" << std::endl;
    LevelBoxData<int, NUMCOMPS, HOST> data_h(layout, Point::Ones(1));
    pout() << "Initializing device data" << std::endl;
    data_d.initialize(d_pointID);
    pout() << "Initializing host data" << std::endl;
    data_h.setVal(-1);

    pout() << "Copying DEVICE --> HOST" << std::endl;
    data_d.copyTo(data_h);
    pout() << "Writing HOST data" << std::endl;
    h5.writeLevel(dx, data_h, "HOST_DATA");
    pout() << "Test complete" << std::endl;
    
    #endif

    #ifdef PR_MPI
    MPI_Finalize();
    #endif
}
