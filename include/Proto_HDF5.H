#pragma once
#include "Proto_LevelBoxData.H"
#include "H5Cpp.h"

namespace Proto
{
    using namespace H5;
    class HDF5Handler
    {
        private:
        
        CompType H5T_PROTO_POINT;
        CompType H5T_PROTO_BOX;
        
        public:
    
        HDF5Handler()
        {
            // define custom types for Proto::Point and Proto::Box
            H5T_PROTO_POINT = CompType(sizeof(Point));
            for (int ii = 0; ii < DIM; ii++)
            {
                char pvar[10];
                sprintf(pvar, "p%i", ii);
                H5T_PROTO_POINT.insertMember(pvar, HOFFSET(Point, m_tuple)+sizeof(int)*ii, PredType::NATIVE_INT);
            }

            H5T_PROTO_BOX = CompType(sizeof(Box));
            H5T_PROTO_BOX.insertMember("low", HOFFSET(Box, m_low), H5T_PROTO_POINT);
            H5T_PROTO_BOX.insertMember("high", HOFFSET(Box, m_high), H5T_PROTO_POINT);
        }
        
        void writeBoxes(std::vector<Box>& a_boxes, const char* a_filename)
        {
            int size = a_boxes.size();
            Box boxes[size];
            for (int ii = 0; ii < size; ii++)
            {
                boxes[ii] = a_boxes[ii];
            }
            H5File file(a_filename, H5F_ACC_TRUNC);
            hsize_t dims[] = {size};
            DataSpace dataspace(1, dims);
            DataSet dataset = file.createDataSet("Boxes", H5T_PROTO_BOX, dataspace);
            dataset.write(boxes, H5T_PROTO_BOX);
        }

    }; // end class HDF5Handler
} // end namespace Proto
