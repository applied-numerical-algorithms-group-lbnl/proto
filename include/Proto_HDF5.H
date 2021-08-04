#pragma once
#include "Proto_LevelBoxData.H"
#include "hdf5.h"

namespace Proto
{
    class HDF5Handler
    {
        public:
    
        HDF5Handler();
        
        template<typename T, unsigned int C, MemType MEM, unsigned int CTR, typename... Args>
        void readLevel(LevelBoxData<T, C, MEM, CTR>& a_data,
                        std::string a_filename,
                        Args... a_params);

        template<typename T, unsigned int C, MemType MEM, unsigned int CTR, typename... Args>
        void writeLevel(const LevelBoxData<T, C, MEM, CTR>& a_data,
                        std::vector<std::string> a_varnames,
                        array<double, DIM>& a_dx,
                        std::string a_filename,
                        Args... a_params);
        
        template<typename T, unsigned int C, MemType MEM, unsigned int CTR, typename... Args>
        void writeLevel(const LevelBoxData<T, C, MEM, CTR>& a_data,
                        std::vector<std::string> a_varnames,
                        std::string a_filename,
                        Args... a_params);
        
        template<typename T, unsigned int C, MemType MEM, unsigned int CTR, typename... Args>
        void writeLevel(const LevelBoxData<T, C, MEM, CTR>& a_data,
                        std::string a_filename,
                        Args... a_params);

        private:
      
        hid_t H5T_PROTO_POINT;
        hid_t H5T_PROTO_BOX; 
        hid_t H5T_PROTO_STRING;

        template<typename T, unsigned int C, MemType MEM, unsigned int CTR>
        void addLevel(hid_t* a_file,
            const LevelBoxData<T, C, MEM, CTR>& a_data,
            array<double, DIM>& a_dx0,
            int a_level);
       /* 
        void writeData(hid_t* a_file,
            std::string a_dataName,
            hid_t* a_space,
            double* a_data);
        
        void writeData(hid_t* a_file,
            std::string a_dataName,
            hid_t* a_space,
            int* a_data);
       */
        /*
        template<typename T>
        void addDataSet(hid_t* a_dataSet,
            hid_t* a_file,
            std::string a_spaceName,
            hid_t* a_dataSpace);
        
        template<typename T>
        void writeData(hid_t* a_file,
            hid_t* a_dataSet,
            hid_t* a_domainSpace,
            hid_t* a_localSpace,
            T* a_data);
        */
        template<typename T>
        void getH5DataType(hid_t* a_type);
        
    }; // end class HDF5Handler
#include "implem/Proto_HDF5Implem.H"
} // end namespace Proto
