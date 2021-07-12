#include "Proto.H"
#include "MHD_Output_Writer.H"
#include "CommonTemplates.H"
#include "Proto_Timer.H"
#include "Proto_WriteBoxData.H"
#include "Proto_LevelBoxData.H"
#include "Proto_ProblemDomain.H"

namespace MHD_Output_Writer {

/**
 * @brief Write out all components of data to filename.vtk. This routine is a no-op for any process
 * that is not process 0.
 * @param data LevelBoxData that is defined on a single box assigned to process 0
 */


	void WriteSinglePatchLevelData(LevelBoxData<double,DIM+NUMCOMPS>& out_data,
	                               const double dx,
	                               const double dy,
	                               const double dz,
	                               const string& filename_data)
	{
		if(procID()==0)
		{
			DataIterator dit=out_data.begin();
#if DIM == 1
			const char* varnames[5];
			varnames[0] = "X";
			varnames[1] = "density";
			varnames[2] = "Vx";
			varnames[3] = "p";
			varnames[4] = "Bx";
#endif

#if DIM == 2
			const char* varnames[8];
			varnames[0] = "X";
			varnames[1] = "Y";
			varnames[2] = "density";
			varnames[3] = "Vx";
			varnames[4] = "Vy";
			varnames[5] = "p";
			varnames[6] = "Bx";
			varnames[7] = "By";
#endif

#if DIM == 3
			const char* varnames[11];
			varnames[0] = "X";
			varnames[1] = "Y";
			varnames[2] = "Z";
			varnames[3] = "density";
			varnames[4] = "Vx";
			varnames[5] = "Vy";
			varnames[6] = "Vz";
			varnames[7] = "p";
			varnames[8] = "Bx";
			varnames[9] = "By";
			varnames[10]= "Bz";
#endif

			double origin[DIM];
			for (int ii = 0; ii < DIM; ii++)
			{
				origin[ii] = 0.0;
			}
			WriteBoxData(filename_data.c_str(),out_data[*dit],varnames,origin,dx);
		}


	}




	void WriteSinglePatchLevelData_nocoord(LevelBoxData<double,NUMCOMPS>& out_data,
	                                       const double dx,
	                                       const double dy,
	                                       const double dz,
	                                       const string& filename_data)
	{
		if(procID()==0)
		{
			DataIterator dit=out_data.begin();
#if DIM == 1
			const char* varnames[4];
			varnames[0] = "density";
			varnames[1] = "Vx";
			varnames[2] = "p";
			varnames[3] = "Bx";
#endif

#if DIM == 2
			const char* varnames[6];
			varnames[0] = "density";
			varnames[1] = "Vx";
			varnames[2] = "Vy";
			varnames[3] = "p";
			varnames[4] = "Bx";
			varnames[5] = "By";
#endif

#if DIM == 3
			const char* varnames[8];
			varnames[0] = "density";
			varnames[1] = "Vx";
			varnames[2] = "Vy";
			varnames[3] = "Vz";
			varnames[4] = "p";
			varnames[5] = "Bx";
			varnames[6] = "By";
			varnames[7]= "Bz";
#endif

			double origin[DIM];
			for (int ii = 0; ii < DIM; ii++)
			{
				origin[ii] = 0.0;
			}
			WriteBoxData(filename_data.c_str(),out_data[*dit],varnames,origin,dx);
		}


	}
}
