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
        void writeLevel(LevelBoxData<T, C, MEM, CTR>& a_data,
                        array<double, DIM>& a_dx,
                        std::string a_filename,
                        Args... a_params);

        private:
      
        hid_t H5T_PROTO_POINT;
        hid_t H5T_PROTO_BOX; 
        hid_t H5T_PROTO_STRING;
            
        template<typename T, unsigned int C, MemType MEM, unsigned int CTR>
        void addLevel(hid_t* a_file,
                      LevelBoxData<T, C, MEM, CTR>& a_data,
                      array<double, DIM>& a_dx0,
                      int a_level);
        
        void writeData(hid_t* a_file,
            std::string a_dataName,
            hid_t* a_space,
            double* a_data);
        
        void writeData(hid_t* a_file,
            std::string a_dataName,
            hid_t* a_space,
            int* a_data);
        
    }; // end class HDF5Handler
#include "implem/Proto_HDF5Implem.H"
} // end namespace Proto
