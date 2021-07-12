#pragma once

#include "Proto.H"
#include "MHDOp.H"

namespace MHD_Output_Writer {
        void WriteSinglePatchLevelData(LevelBoxData<double,DIM+NUMCOMPS>& out_data,
                                       const double dx,
                                       const double dy,
                                       const double dz,
                                       const string& filename_data);

        void WriteSinglePatchLevelData_nocoord(LevelBoxData<double,NUMCOMPS>& out_data,
                                               const double dx,
                                               const double dy,
                                               const double dz,
                                               const string& filename_data);
}
