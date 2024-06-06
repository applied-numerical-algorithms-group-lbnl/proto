#include "Proto.H"
#include "MHDReader.H"
#define NGHOST 6
using namespace Proto;

namespace {
    template<typename T, unsigned int C, MemType MEM>
    MBAMRGrid regrid(MBAMRData<T,C,MEM>& data, T threshold)
    {
        MBAMRGrid grid = data.grid();
        int finestLevel = grid.numLevels()-1;
        std::vector<MBPatchID_t> patches;
        auto finestGrid = grid.getLevel(finestLevel);
        for (auto iter : finestGrid)
        {
            auto& di = data.getLevel(finestLevel)[iter];
            BoxData<T,C,MEM> tmp(finestGrid[iter]);
            di.copyTo(tmp);
            T value = tmp.absMax(0);
            if (value > threshold)
            {
                patches.push_back(MBPatchID_t(finestGrid.point(iter), finestGrid.block(iter)));
            }
        }
        grid.setPatches(finestLevel, patches);
        for (int lvl = finestLevel-1; lvl > 0; lvl--)
        {
            for (BlockIndex bi = 0; bi < grid.numBlocks(); bi++)
            {
                grid.getBlock(bi).enforceNesting(lvl);
            }
        }
        return grid;
    }

}

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
    Point ghost = Point::Ones(NGHOST);
    ghost[rCoord] = 1;
    
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
    int domainSize = 64;
    int boxSize = 16;
    int thickness = 1;
    
    auto domain = CubedSphereShell::Domain(domainSize, thickness, rCoord);
    Point boxSizeVect = Point::Ones(boxSize);
    boxSizeVect[rCoord] = thickness;
    MBDisjointBoxLayout layout(domain, boxSizeVect);
    MBLevelBoxData<double, 8, HOST> dstData(layout, ghost);
    //MBLevelBoxData<double, 8, HOST> dstData2(layout, Point::Basis(rCoord) + NGHOST*Point::Basis(thetaCoord) + NGHOST*Point::Basis(phiCoord));
    MBLevelBoxData<double, 8, HOST> dstData2(layout, ghost);
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

    // used to generate a figure
    // Gaussian filter
#if 0
    Stencil<double> S = 41.0*Shift::Zeros();
    Box B1 = Box::Kernel(1).grow(rCoord, -1);
    Box B2 = Box::Kernel(2).grow(rCoord, -2);
    for (auto si : B2)
    {
        if (B1.contains(si)) { continue; }
        switch (si.abs().sum())
        {
            case 2: S += 7.0*Shift(si); break;
            case 3: S += 4.0*Shift(si); break;
            case 4: S += 1.0*Shift(si); break;
        }
    }
    for (auto si : B1)
    {
        if (si == Point::Ones()) { continue; }
        switch (si.abs().sum())
        {
            case 1: S += 26.0*Shift(si);
            case 2: S += 16.0*Shift(si);
        }
    }
    S *= (1.0/241.0);

    Point refRatio = Point::Ones(4);
    refRatio[rCoord] = 1;
    int numLevels = 2;
    MBAMRGrid grid(domain, boxSizeVect, refRatio, numLevels);
    MBAMRData<double, 8, HOST> amrData(grid, ghost);
    MBAMRData<double, 8, HOST> filteredAMRData(grid, Point::Zeros());
    MBAMRMap<MBMap_CubedSphereShell, HOST> amrMap(grid, ghost);
    filteredAMRData.setVal(-1);
    MBLevelBoxData<double, 8, HOST> filteredData(layout, Point::Zeros());

    for (auto iter : layout)
    {
        auto& fi = filteredData[iter];
        auto& si = dstData2[iter];
        fi |= S(si);
    }
    //dstData2.copyTo(amrData.getLevel(0));
    filteredData.copyTo(amrData.getLevel(0));

    amrData.interpolate();

    for (int li = 0; li < numLevels; li++)
    {
        for (auto iter : grid.getLevel(li))
        {
            auto& fi = filteredAMRData.getLevel(li)[iter];
            auto& si = amrData.getLevel(li)[iter];
            fi |= S(si);
        }
    }
    double threshold = 500;
    auto newGrid = regrid(amrData, threshold);

    MBAMRData<double, 8, HOST> newAMRData(newGrid, ghost);
    MBAMRMap<MBMap_CubedSphereShell, HOST> newAMRMap(newGrid, ghost);
    newAMRData.setVal(7);

    h5.writeMBLevel(map, filteredData, "FILTERED_DATA");
    h5.writeMBAMRData({"1","2","3","4","5","6","7","8"}, newAMRMap, newAMRData, "REGRID_AMR_DATA");
    h5.writeMBAMRData({"1","2","3","4","5","6","7","8"}, amrMap, filteredAMRData, "FILTERED_AMR_DATA");
    h5.writeMBAMRData({"density","2","3","4","5","6","7","8"}, amrMap, amrData, "MHD_STATE");
#endif
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
