#include "Proto.H"
#include "hdf5.h"

using namespace Proto;

class MHDReader
{
    public:
    inline MHDReader(){}
    
    template<typename T, MemType MEM>
    inline void file_to_BC(MBLevelBoxData<T,8,MEM>& a_dstData, 
                           const MBLevelMap<MBMap_CubedSphereShell, MEM>& a_map,
                           std::string a_filename,
                           T a_time)
    {
        int rCoord = CUBED_SPHERE_SHELL_RADIAL_COORD;
        int thetaCoord = (rCoord + 1) % 3;
        int phiCoord = (rCoord + 2) % 3;
        
        // Read meta data from file
        std::vector<double> dtheta;
        Box B = readDomain(a_filename);
        readGeom(dtheta, a_filename);
        BoxData<double, 8, HOST> data_0(B);
        readData(data_0, a_filename);
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
        data_0.rotate(BT,R);

        // Extend the box to double the size in the phi direction
        Box bx1 = data_0.box();
        Point high = bx1.high();
        int domainSizez = high[phiCoord]+1;
        high[phiCoord] = 2*(bx1.high()[phiCoord]+1)-1;
        Point low = bx1.low();
        low[phiCoord] = -1;
        Box bx2(low, high);
        BoxData<double, 8, HOST> data_1(bx2);
        BoxData<double, 8, HOST> data_2(bx2);
        data_0.copyTo(data_1,bx1,Point::Basis(phiCoord)*(bx1.high()[phiCoord]+1));
        data_0.copyTo(data_1);

        // Rotate the data
        double carr_rot_time = 25.38*24*60*60; // Seconds
        // We should use sidereal time for this. 25.38 days. That's the rotation time from a fixed location.
        // Carrington rotation time (27.2753 days) is from Earth's prespective.
        double angle_to_rotate = fmod(360*a_time/carr_rot_time,360);
        int cells_to_rotate = angle_to_rotate/(360/domainSizez);
        double needed_fraction = angle_to_rotate/(360/domainSizez) - cells_to_rotate;
        cells_to_rotate = cells_to_rotate % domainSizez;
        cells_to_rotate = domainSizez - cells_to_rotate;
        static Stencil<double> m_right_shift;
        m_right_shift = (1.0-needed_fraction)*Shift(Point::Zeros()) + (needed_fraction)*Shift(-Point::Basis(phiCoord));
        data_1.copyTo(data_2,bx2,Point::Basis(phiCoord)*(-cells_to_rotate));
        data_1 = m_right_shift(data_2);
        data_1.copyTo(data_0,bx1,Point::Zeros());

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

        // MBLevelBoxData<double, 8, HOST> dstData(a_JU.layout(), Point::Basis(rCoord) + NGHOST*Point::Basis(thetaCoord) + NGHOST*Point::Basis(phiCoord));
        // a_dstData.define(a_layout, {Point::Zeros(),Point::Zeros(),Point::Zeros(),Point::Zeros()});
        a_dstData.setVal(0);
        
        // Get the interpolation or Nearest operator
        auto BCOp = CubedSphereShell::BCNearestOp(a_map, srcLayout, dtheta, Side::Lo, NGHOST);

        // Interpolate data
        BCOp.apply(a_dstData, data);
    }
     
    template<typename... Args>
    inline void readFile(hid_t* a_file, std::string a_filename, Args... a_params)
    {
        char fname[100];
        // if (a_filename.substr(a_filename.find_last_of(".") + 1) != "hdf5")
        // {
        //     a_filename += ".hdf5";
        // }

        // The following pragmas suppress some GCC warnings that are overly conservative
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wformat-security"
        sprintf(fname, a_filename.c_str(), a_params...);
#pragma GCC diagnostic pop

#ifdef PR_MPI
        auto mpi_comm = MPI_COMM_WORLD;
        auto mpi_info = MPI_INFO_NULL;

        auto p_access = H5Pcreate(H5P_FILE_ACCESS);
        assert(H5Pset_fapl_mpio(p_access, mpi_comm, mpi_info) >= 0);
        *a_file = H5Fopen(fname, H5F_ACC_RDONLY, p_access);
        assert(H5Pclose(p_access) >= 0);
#else
        *a_file = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
#endif
    }
   
    inline Box _readDomain(hid_t* a_file)
    {
        // READ NUM COMPONENTS
        //int numComps;
        //auto att_comps = H5Aopen(*a_file, "num_components", H5P_DEFAULT);
        //assert(H5Aread(att_comps, H5T_STD_I32LE, &numComps) >= 0); 
        //assert(H5Aclose(att_comps) >= 0);

        //PROTO_ASSERT(numComps == C,
        //        "MHDReader::readLevel | Error: Incorrect number of components. %i != %i",
        //        C, numComps);

        // READ DOMAIN SIZES
        int domainSize[2];
        auto att_domain = H5Aopen(*a_file, "domain", H5P_DEFAULT);
        assert(H5Aread(att_domain, H5T_STD_I32LE, domainSize) >= 0); 
        assert(H5Aclose(att_domain) >= 0);
        
        int rCoord = CUBED_SPHERE_SHELL_RADIAL_COORD;
        int thetaCoord = (rCoord + 1) % 3;
        int phiCoord = (rCoord + 2) % 3;
        
        // these values are backwards in the data file
        Point sizes = Point::Ones();
        sizes[thetaCoord] = domainSize[0];
        sizes[phiCoord] = domainSize[1];

        return Box(sizes);
    }
    
    template<typename... Args>
    inline Box readDomain(
            std::string a_filename,
            Args... a_params)
    {
        hid_t file;
        readFile(&file, a_filename, a_params...);
        Box domainBox = _readDomain(&file);
        barrier();
        assert(H5Fclose(file) >= 0);
        return domainBox;
    }
    
    template<typename T, unsigned int C, MemType MEM, typename... Args>
    inline void readDataset(hid_t* a_file,
            BoxData<T, C, MEM>& a_data, 
            int a_index)
    {
        hid_t H5T_T;
        HDF5Handler::getH5DataType<T>(&H5T_T);
        auto ds_data = H5Dopen2(*a_file, ("/data" + std::to_string(a_index)).c_str(), H5P_DEFAULT);
        if (MEM == HOST)
        {
            assert(H5Dread(ds_data, H5T_T, H5S_ALL, H5S_ALL, H5P_DEFAULT, a_data.data()) >= 0);
        } else {
            int dataSize = a_data.size();
            T* rawData = (T*)proto_malloc<HOST>(dataSize);
            assert(H5Dread(ds_data, H5T_T, H5S_ALL, H5S_ALL, H5P_DEFAULT, rawData) >= 0);
            proto_memcpy<HOST, MEM>(rawData, a_data.data(), dataSize);
            proto_free<HOST>(rawData);
        }
        assert(H5Dclose(ds_data) >= 0);
    }
    
    template<typename T, unsigned int C, MemType MEM, typename... Args>
    inline void readData(
            BoxData<T, C, MEM>& a_data,
            std::string a_filename,
            Args... a_params)
    {
        hid_t file;
        readFile(&file, a_filename, a_params...);
        Box domainBox = _readDomain(&file);
        a_data.define(domainBox);
        readDataset(&file, a_data, 0); 
        barrier();
        assert(H5Fclose(file) >= 0);
    }

    template<typename T, unsigned int C, MemType MEM, typename... Args>
    inline void readData(
            std::vector<BoxData<T, C, MEM>>& a_data,
            std::string a_filename,
            Args... a_params)
    {
        hid_t file;
        readFile(&file, a_filename, a_params...);
        Box domainBox = _readDomain(&file);

        int numDatasets;
        auto att_domain = H5Aopen(file, "num_datasets", H5P_DEFAULT);
        assert(H5Aread(att_domain, H5T_STD_I32LE, &numDatasets) >= 0); 
        assert(H5Aclose(att_domain) >= 0);
        
        a_data.clear();
        a_data.resize(numDatasets);
        for (int ii = 0; ii < numDatasets; ii++)
        {
            a_data[ii].define(domainBox);
            readDataset(&file, a_data[ii], ii);
        }
        barrier();
        assert(H5Fclose(file) >= 0);
    }
    
    template<typename T, typename... Args>
    inline void readGeom(
            std::vector<T>& a_data,
            std::string a_filename,
            Args... a_params)
    {
        PROTO_ASSERT(DIM == 3,
                "MHDReader::readGeom | Error: Invalid for DIM != 3");

        PR_TIME("MHDReader::readGeom");
        hid_t file;
        readFile(&file, a_filename, a_params...);
        
        hid_t H5T_T;
        HDF5Handler::getH5DataType<T>(&H5T_T);

        // READ DOMAIN SIZES
        int domainSize[2];
        auto att_domain = H5Aopen(file, "domain", H5P_DEFAULT);
        assert(H5Aread(att_domain, H5T_STD_I32LE, domainSize) >= 0); 
        assert(H5Aclose(att_domain) >= 0);
       
        // READ DATA
        int dataSize = domainSize[1];
        auto ds_data = H5Dopen2(file, "/geometry/dtheta", H5P_DEFAULT);

        // DEFINE VECTOR
        a_data.resize(dataSize);
        
        // READ DATA
        assert(H5Dread(ds_data, H5T_T, H5S_ALL, H5S_ALL, H5P_DEFAULT, &a_data[0]) >= 0);
        
        // CLEAN UP
        assert(H5Dclose(ds_data) >= 0);
        barrier();
        assert(H5Fclose(file) >= 0);
    }
};
