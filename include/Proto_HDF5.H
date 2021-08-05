#pragma once
#ifndef PROTO_HDF5_H
#define PROTO_HDF5_H
#include "Proto_LevelBoxData.H"
#include "hdf5.h"

namespace Proto
{
    class HDF5Handler
    {
        public:
    
        inline HDF5Handler();
        
        template<typename T, unsigned int C, MemType MEM, unsigned int CTR, typename... Args>
        inline void readLevel(LevelBoxData<T, C, MEM, CTR>& a_data,
                        std::string a_filename,
                        Args... a_params);

        template<typename T, unsigned int C, MemType MEM, unsigned int CTR, typename... Args>
        inline void writeLevel(const LevelBoxData<T, C, MEM, CTR>& a_data,
                        std::vector<std::string> a_varnames,
                        array<double, DIM>& a_dx,
                        std::string a_filename,
                        Args... a_params);
        
        template<typename T, unsigned int C, MemType MEM, unsigned int CTR, typename... Args>
        inline void writeLevel(const LevelBoxData<T, C, MEM, CTR>& a_data,
                        std::vector<std::string> a_varnames,
                        std::string a_filename,
                        Args... a_params);
        
        template<typename T, unsigned int C, MemType MEM, unsigned int CTR, typename... Args>
        inline void writeLevel(const LevelBoxData<T, C, MEM, CTR>& a_data,
                        std::string a_filename,
                        Args... a_params);

        private:
      
        hid_t H5T_PROTO_POINT;
        hid_t H5T_PROTO_BOX; 
        hid_t H5T_PROTO_STRING;

        template<typename T, unsigned int C, MemType MEM, unsigned int CTR>
        inline void addLevel(hid_t* a_file,
            const LevelBoxData<T, C, MEM, CTR>& a_data,
            array<double, DIM>& a_dx0,
            int a_level);
        
        template<typename T>
        inline void getH5DataType(hid_t* a_type) {}
        
    }; // end class HDF5Handler
#include "implem/Proto_HDF5Implem.H"
} // end namespace Proto
#endif //end of include guard
