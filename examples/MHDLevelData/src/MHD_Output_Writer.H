#pragma once

#include "Proto.H"
#include "MHDOp.H"

namespace MHD_Output_Writer {
        void WriteSinglePatchLevelData(LevelBoxData<double,DIM+NUMCOMPS>& out_data,
	                               const double dx,
	                               const double dy,
	                               const double dz,
                                       const int k,
                                       const double time,
	                               const string& filename_data);

        void WriteSinglePatchLevelData_nocoord(LevelBoxData<double,NUMCOMPS>& out_data,
                                               const double dx,
                                               const double dy,
                                               const double dz,
                                               const string& filename_data);

        void WriteBoxData_array(BoxData<double,DIM+NUMCOMPS>& out_data,
   	                               const double dx,
   	                               const double dy,
   	                               const double dz,
   	                               const string& filename_data);

        void WriteBoxData_array_nocoord(const BoxData<double,NUMCOMPS>& out_data,
   	                                       const double dx,
   	                                       const double dy,
   	                                       const double dz,
   	                                       const string& filename_data);

        void WriteBoxData_Nsd_nocoord(BoxData<double,DIM*DIM>& out_data,
	                                       const double dx,
	                                       const double dy,
	                                       const double dz,
	                                       const string& filename_data);
}
