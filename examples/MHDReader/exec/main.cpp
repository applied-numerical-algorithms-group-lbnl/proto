#include "Proto.H"
#include "MHDReader.H"
#define NGHOST 6
using namespace Proto;

int main(int argc, char** argv)
{
    #ifdef PR_MPI
    MPI_Init(&argc, &argv);
    #endif
    
    MHDReader reader;
    HDF5Handler h5;

#ifdef PR_MMB
#if DIM == 3
    int rCoord = CUBED_SPHERE_SHELL_RADIAL_COORD;
    int thetaCoord = (rCoord + 1) % 3;
    int phiCoord = (rCoord + 2) % 3;
    
    // Read meta data from "DATA.hdf5"
    std::vector<double> dtheta;
    Box B = reader.readDomain("DATA");
    reader.readGeom(dtheta, "DATA");
    BoxData<double, 8, HOST> data_0(B);

    reader.readData(data_0, "DATA");

    MPI_Bcast(data_0.data(), B.size()*8 ,MPI_DOUBLE,0,MPI_COMM_WORLD);

    // Transpose Phi and Theta
    Point LT = B.low(); 
    LT[thetaCoord] = B.low()[phiCoord];
    LT[phiCoord] = B.low()[thetaCoord];
    Point HT = B.high(); 
    HT[thetaCoord] = B.high()[phiCoord];
    HT[phiCoord] = B.high()[thetaCoord];
    Box BT(LT,HT);
    auto t = Point::Basis(thetaCoord);
    auto p = Point::Basis(phiCoord);
    CoordPermutation R{{t,p},{p,t}};
    // h5.writePatch(1, data_0, "SRC_DATA");
    data_0.rotate(BT,R);
    // h5.writePatch(1, data_0, "SRC_DATA_T");

    // Create a single block, single patch MBLevelBoxData
    // You may have to do this on each processor
    MBProblemDomain srcDomain(numProc());
    for (BlockIndex bi = 0; bi < numProc(); bi++)
    {
        srcDomain.defineDomain(bi,BT.sizes());
    }
    MBDisjointBoxLayout srcLayout(srcDomain, BT.sizes());
    MBLevelBoxData<double, 8, HOST> data(srcLayout, Point::Zeros());
    data_0.copyTo(data[*srcLayout.begin()]);

    // Cubed Sphere data
    int domainSize = 128;
    int boxSize = 128;
    int thickness = 1;
    
    auto domain = CubedSphereShell::Domain(domainSize, thickness, rCoord);
    Point boxSizeVect = Point::Ones(boxSize);
    boxSizeVect[rCoord] = thickness;
    MBDisjointBoxLayout layout(domain, boxSizeVect);
    MBLevelBoxData<double, 8, HOST> dstData(layout, Point::Basis(rCoord) + NGHOST*Point::Basis(thetaCoord) + NGHOST*Point::Basis(phiCoord));
    MBLevelBoxData<double, 8, HOST> dstData2(layout, Point::Basis(rCoord) + NGHOST*Point::Basis(thetaCoord) + NGHOST*Point::Basis(phiCoord));
    dstData.setVal(0);

    auto map = CubedSphereShell::Map(layout, Point::Ones());
    
    // Get the interpolation operator
    // auto interpOp = CubedSphereShell::BCInterpOp(map, srcLayout, dtheta, Side::Lo);
    auto interpOp = CubedSphereShell::BCNearestOp(map, srcLayout, dtheta, Side::Lo, NGHOST);

    // Interpolate data
    interpOp.apply(dstData, data);

    for (auto dit : dstData.layout())
    {
    Point source_lo = Point(-1,dstData[dit].box().low()[1],dstData[dit].box().low()[2]);
	Point source_hi = Point(-1,dstData[dit].box().high()[1],dstData[dit].box().high()[2]);
    Box sourceBox(source_lo,source_hi);
    dstData[dit].copyTo(dstData2[dit],sourceBox,Point::Basis(0)*(1));
    }

    h5.writeMBLevel(data, "SRC_DATA");
    h5.writeMBLevel(map, dstData, "DST_DATA");
    h5.writeMBLevel(map, dstData2, "DST_DATA2");
#else
    std::cout << "No test run. DIM must be equal to 3" << std::endl;
#endif
#else
    // Read data from "DATA.hdf5"
    BoxData<double, 8, HOST> data;
    std::vector<BoxData<double, 8, HOST>> timeData;
    std::vector<double> dtheta;
    reader.readData(data, "DATA");
    reader.readData(timeData, "TIME_DATA");
    reader.readGeom(dtheta, "DATA");

    // Build grids
    Box domainBox(Point(180,360,180));
    Point boxSize = Point::Ones(90);
    std::array<bool, DIM> periodicity = {false, true, false}; //not read from HDF5 file
    ProblemDomain domain(domainBox, periodicity);
    DisjointBoxLayout layout(domain, boxSize);

    // Initialize data
    LevelBoxData<double, 8, HOST> state(layout, Point::Ones());
    state.setToZero(); // replace with initialization

    // Apply boundary condition
    for (auto index : layout)
    {
        auto& state_i = state[index];
        Point s(-1,0,0);            // shift in -x direction to put data in ghost region
        Box box_i = state_i.box().shift(s) & data.box(); //region of data to copy
        data.copyTo(state_i, box_i, s); 
    }

    h5.writeLevel(1, state, "STATE");

    for (int ii = 0; ii < timeData.size(); ii++)
    {
        h5.writePatch(1, timeData[ii], "STATE_%i", ii);
    }
    //do something with dtheta...
#endif
    #ifdef PR_MPI
    MPI_Finalize();
    #endif
}
