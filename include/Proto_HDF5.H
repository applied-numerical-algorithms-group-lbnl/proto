#pragma once
#include "Proto_LevelBoxData.H"
#include "H5Cpp.h"

namespace Proto
{
    using namespace H5;
    class HDF5Handler
    {
        public:
    
        HDF5Handler();
        void writeBoxes(std::vector<Box>& a_boxes, const char* a_filename);

        template<typename T, unsigned int C, MemType MEM, unsigned int CTR, typename... Args>
        void writeLevel(LevelBoxData<T, C, MEM, CTR>& a_data,
                        array<double, DIM>& a_dx,
                        std::string a_filename,
                        Args... a_params);

        private:
        
        CompType H5T_PROTO_POINT;
        CompType H5T_PROTO_BOX;
        StrType  H5T_PROTO_STRING;
                                
        template<typename T, unsigned int C, MemType MEM, unsigned int CTR>
        void addLevel(H5File& a_file,
                      LevelBoxData<T, C, MEM, CTR>& a_data,
                      array<double, DIM>& a_dx0,
                      int a_level);
        
        void writeData(H5File& a_file,
            std::string a_dataName,
            DataSpace& a_space,
            double* a_data);
        
        void writeData(H5File& a_file,
            std::string a_dataName,
            DataSpace& a_space,
            int* a_data);
        
    }; // end class HDF5Handler
#include "implem/Proto_HDF5Implem.H"
} // end namespace Proto
