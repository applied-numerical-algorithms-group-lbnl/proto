#include "Proto.H"
#include "hdf5.h"

using namespace Proto;

class MHDReader
{
    public:
    inline MHDReader(){}

    template<typename T, unsigned int C, MemType MEM, typename... Args>
    inline void readData(
            BoxData<T, C, MEM>& a_data,
            std::string a_filename,
            Args... a_params)
    {
        PROTO_ASSERT(DIM == 3,
                "MHDReader::readLevel | Error: Invalid for DIM != 3");

        PR_TIME("MHDReader::readLevel");
        // PARSE THE FILE NAME
        char fname[100];
        if (a_filename.substr(a_filename.find_last_of(".") + 1) != "hdf5")
        {
            a_filename += ".hdf5";
        }

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
        auto file = H5Fopen(fname, H5F_ACC_RDONLY, p_access);
        assert(H5Pclose(p_access) >= 0);
#else
        auto file = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
#endif
        hid_t H5T_T;
        HDF5Handler::getH5DataType<T>(&H5T_T);

        // READ DOMAIN SIZES
        int domainSize[2];
        auto att_domain = H5Aopen(file, "domain", H5P_DEFAULT);
        assert(H5Aread(att_domain, H5T_STD_I32LE, domainSize) >= 0); 
        assert(H5Aclose(att_domain) >= 0);
        
        // READ NUM COMPONENTS
        int numComps;
        auto att_comps = H5Aopen(file, "num_components", H5P_DEFAULT);
        assert(H5Aread(att_comps, H5T_STD_I32LE, &numComps) >= 0); 
        assert(H5Aclose(att_comps) >= 0);

        PROTO_ASSERT(numComps == C,
                "MHDReader::readLevel | Error: Incorrect number of components. %i != %i",
                C, numComps);

        // DEFINE BOX DATA
        Box domainBox(Point(1, domainSize[0], domainSize[1]));
        a_data.define(domainBox);
        
        // READ DATA
        int dataSize = domainSize[0]*domainSize[1]*numComps;
        auto ds_data = H5Dopen2(file, "/data0", H5P_DEFAULT);
        if (MEM == HOST)
        {
            assert(H5Dread(ds_data, H5T_T, H5S_ALL, H5S_ALL, H5P_DEFAULT, a_data.data()) >= 0);
        } else {
            T* rawData = (T*)proto_malloc<HOST>(dataSize);
            assert(H5Dread(ds_data, H5T_T, H5S_ALL, H5S_ALL, H5P_DEFAULT, rawData) >= 0);
            proto_memcpy<HOST, MEM>(rawData, a_data.data(), dataSize);
            proto_free<HOST>(rawData);
        }
        assert(H5Dclose(ds_data) >= 0);
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
        // PARSE THE FILE NAME
        char fname[100];
        if (a_filename.substr(a_filename.find_last_of(".") + 1) != "hdf5")
        {
            a_filename += ".hdf5";
        }

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
        auto file = H5Fopen(fname, H5F_ACC_RDONLY, p_access);
        assert(H5Pclose(p_access) >= 0);
#else
        auto file = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
#endif
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
