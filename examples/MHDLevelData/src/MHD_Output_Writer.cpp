#include "Proto.H"
#include "MHD_Output_Writer.H"
#include "CommonTemplates.H"
#include "Proto_Timer.H"
#include "Proto_WriteBoxData.H"
#include "Proto_LevelBoxData.H"
#include "Proto_ProblemDomain.H"
#include "MHD_Mapping.H"
#include "MHD_Input_Parsing.H"
extern Parsefrominputs inputs;

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
								   const int k,
								   const double time,
	                               const string& filename_data)
	{
		if(procID()==0)
		{
			DataIterator dit=out_data.begin();

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
			// WriteBoxData(filename_data.c_str(),out_data[*dit],varnames,origin,1);
		}
		array<double, DIM> a_dx;
		a_dx[0] = dx;
		a_dx[1] = dy;
		a_dx[2] = dz;
		HDF5Handler h5;

#if DIM == 2		
		h5.writeLevel({"X","Y","density","Vx","Vy", "p","Bx","By"}, 1, out_data, filename_data);
#endif
#if DIM == 3	
		double step = 1.0*k;
		h5.setTime(time);	
		h5.setTimestep(step);	
		h5.writeLevel({"X","Y","Z","density","Vx","Vy","Vz", "p","Bx","By","Bz"}, 1, out_data, filename_data);
#endif


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
			// WriteBoxData(filename_data.c_str(),out_data[*dit],varnames,origin,1);
			
		}
		HDF5Handler h5;
#if DIM == 2		
		h5.writeLevel({"density","Vx","Vy", "p","Bx","By"}, 1, out_data, filename_data);
#endif
#if DIM == 3		
		h5.writeLevel({"density","Vx","Vy","Vz", "p","Bx","By","Bz"}, 1, out_data, filename_data);
#endif
	}




	void WriteBoxData_array(BoxData<double,DIM+NUMCOMPS>& out_data,
	                               const double dx,
	                               const double dy,
	                               const double dz,
	                               const string& filename_data)
	{
		if(procID()==0)
		{

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
			WriteBoxData(filename_data.c_str(),out_data,varnames,origin,1);
		}

	}

	void WriteBoxData_array_nocoord(const BoxData<double,NUMCOMPS>& out_data,
	                                       const double dx,
	                                       const double dy,
	                                       const double dz,
	                                       const string& filename_data)
	{
		Box dbx0 = out_data.box();
		BoxData<double,DIM> eta(dbx0), X_center(dbx0);
		MHD_Mapping::eta_calc(eta,dbx0,dx, dy, dz);
		// MHD_Mapping::eta_to_x_calc(X_center,eta);
		if (inputs.grid_type_global == 2){
			// MHD_Mapping::eta_to_x_ave_calc(X_center,eta,dx,dy,dz);  // This is not working for case 16
			MHD_Mapping::eta_to_x_calc(X_center,eta);
		} else { 
			MHD_Mapping::eta_to_x_calc(X_center,eta);
		}
		BoxData<double,NUMCOMPS+DIM> out_data2;
		MHD_Mapping::out_data_calc(out_data2,X_center,out_data);

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
			WriteBoxData(filename_data.c_str(),out_data2,varnames,origin,1);
	}


	void WriteBoxData_Nsd_nocoord(BoxData<double,DIM*DIM>& out_data,
	                                       const double dx,
	                                       const double dy,
	                                       const double dz,
	                                       const string& filename_data)
	{
		Box dbx0 = out_data.box();
		BoxData<double,DIM> eta(dbx0), X_center(dbx0);
		MHD_Mapping::etaFace_calc(eta,dbx0,dx, dy, dz,0);
		MHD_Mapping::eta_to_x_calc(X_center,eta);
		BoxData<double,DIM*DIM+DIM> out_data2;
		MHD_Mapping::out_data2_calc(out_data2,X_center,out_data);

#if DIM == 2
			const char* varnames[6];
			varnames[0] = "X";
			varnames[1] = "Y";
			varnames[2] = "N11";
			varnames[3] = "N21";
			varnames[4] = "N12";
			varnames[5] = "N22";
#endif

#if DIM == 3
			const char* varnames[12];
			varnames[0] = "X";
			varnames[1] = "Y";
			varnames[2] = "Z";
			varnames[3] = "N11";
			varnames[4] = "N21";
			varnames[5] = "N31";
			varnames[6] = "N12";
			varnames[7] = "N22";
			varnames[8] = "N32";
			varnames[9] = "N13";
			varnames[10]= "N23";
			varnames[11]= "N33";

#endif

			double origin[DIM];
			for (int ii = 0; ii < DIM; ii++)
			{
				origin[ii] = 0.0;
			}
			WriteBoxData(filename_data.c_str(),out_data2,varnames,origin,1);
	}



}
