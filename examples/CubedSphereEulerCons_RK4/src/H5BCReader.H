#include "Proto.H"
#include "hdf5.h"
// #include "H5BC_Input_Parsing.H"
// extern Parsefrominputs inputs;

using namespace Proto;

class H5BCReader
{
    public:
    inline H5BCReader(){}
   
     
    template<typename... Args>
    inline void readFile(hid_t* a_file, std::string a_filename, Args... a_params)
    {
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
        *a_file = H5Fopen(fname, H5F_ACC_RDONLY, p_access);
        assert(H5Pclose(p_access) >= 0);
#else
        *a_file = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
#endif
    }
   
    template<unsigned int C>
    inline Box readDomain(hid_t* a_file)
    {
        // READ NUM COMPONENTS
        int numComps;
        auto att_comps = H5Aopen(*a_file, "num_components", H5P_DEFAULT);
        assert(H5Aread(att_comps, H5T_STD_I32LE, &numComps) >= 0); 
        assert(H5Aclose(att_comps) >= 0);

        PROTO_ASSERT(numComps == C,
                "H5BCReader::readLevel | Error: Incorrect number of components. %i != %i",
                C, numComps);

        // READ DOMAIN SIZES
        int domainSize[2];
        auto att_domain = H5Aopen(*a_file, "domain", H5P_DEFAULT);
        assert(H5Aread(att_domain, H5T_STD_I32LE, domainSize) >= 0); 
        assert(H5Aclose(att_domain) >= 0);
        
        return Box(Point(1, domainSize[0], domainSize[1]));
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
            std::vector<BoxData<T, C, MEM>>& a_data,
            std::string a_filename,
            Args... a_params)
    {
        PROTO_ASSERT(DIM == 3,
                "H5BCReader::readLevel | Error: Invalid for DIM != 3");

        PR_TIME("H5BCReader::readLevel");
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
                "H5BCReader::readLevel | Error: Incorrect number of components. %i != %i",
                C, numComps);

        // int read_every = domainSize[0]/inputs.domainSizez;
        int read_every = 1;

        int numDatasets;
        att_domain = H5Aopen(file, "num_datasets", H5P_DEFAULT);
        assert(H5Aread(att_domain, H5T_STD_I32LE, &numDatasets) >= 0); 
        assert(H5Aclose(att_domain) >= 0);
        
        // READ DATA
        int dataSize = domainSize[0]*domainSize[1]*numComps;
        
        auto ds_data = H5Dopen2(file, "/data0", H5P_DEFAULT);

        a_data.clear();
        a_data.resize(numDatasets);
        for (int ii = 0; ii < numDatasets; ii++)
        {
            ds_data = H5Dopen2(file, ("/data" + std::to_string(ii)).c_str(), H5P_DEFAULT);
            // DEFINE BOX DATA
            Box domainBox(Point::Zeros(),Point(0, domainSize[1]/read_every-1, 2*domainSize[0]/read_every-1));
            a_data[ii].define(domainBox);
            
            Box domainBox2(Point::Zeros(),Point(0, domainSize[0]-1, domainSize[1]-1));
            BoxData<double,C> a_data2(domainBox2);
            
            
            if (MEM == HOST)
            {
            //     assert(H5Dread(ds_data, H5T_T, H5S_ALL, H5S_ALL, H5P_DEFAULT, a_data.data()) >= 0);
                assert(H5Dread(ds_data, H5T_T, H5S_ALL, H5S_ALL, H5P_DEFAULT, a_data2.data()) >= 0);
                // This is needed to make transpose. h5 file is (phi, theta). We need (theta, phi).
                for (int c = 0; c < numComps; c++){
                    for (int i = 0; i < domainSize[0]/read_every; i++ )
                    {
                        for (int j = 0; j < domainSize[1]/read_every; j++ )
                        {
                            // Index in the final matrix.
                            int index1 = 2*c*domainSize[0]/read_every*domainSize[1]/read_every + i*domainSize[1]/read_every+j;
                            // Index in the original matrix.
                            int index2 = c*domainSize[0]*domainSize[1] + (j*domainSize[0]+i)*read_every;
                            a_data[ii].data()[index1] = a_data2.data()[index2];
                            a_data[ii].data()[index1+domainSize[0]/read_every*domainSize[1]/read_every] = a_data[ii].data()[index1];
                        }
                    }   
                }

            } else {
                T* rawData = (T*)proto_malloc<HOST>(dataSize);
                assert(H5Dread(ds_data, H5T_T, H5S_ALL, H5S_ALL, H5P_DEFAULT, rawData) >= 0);
                proto_memcpy<HOST, MEM>(rawData, a_data[ii].data(), dataSize);
                proto_free<HOST>(rawData);
            }
        }
        assert(H5Dclose(ds_data) >= 0);
        barrier();
        // assert(H5Fclose(file) >= 0);
    }


    
    template<typename T, typename... Args>
    inline void readGeom(
            std::vector<T>& a_data,
            std::string a_filename,
            Args... a_params)
    {
        PROTO_ASSERT(DIM == 3,
                "H5BCReader::readGeom | Error: Invalid for DIM != 3");

        PR_TIME("H5BCReader::readGeom");
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
