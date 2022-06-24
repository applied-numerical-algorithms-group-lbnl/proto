#include <cstdio>
#include <cstring>
#include <cassert>
#include <cmath>

#include <vector>
#include <memory>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "Proto.H"
#include "MHDLevelDataRK4.H"
#include "Proto_WriteBoxData.H"
#include "Proto_Timer.H"

using namespace std;
using namespace Proto;


double time_globalll;
int grid_type_global;
int LowBoundType;
int HighBoundType;
bool limiter_apply;
bool slope_flattening_apply;
bool non_linear_visc_apply;
bool linear_visc_apply;
bool pole_correction;
int init_condition_type;
int Riemann_solver_type;
int domainSizex;
int domainSizey;
int domainSizez;
int BoxSize;
int lim_count;
int lim_count2;
double current_iter;
bool initialize_in_spherical_coords;
// void InitializeMHDLevelDataState(MHDLevelDataState& state)
// {
	
// }

typedef BoxData<double,1,HOST> Scalar;
typedef BoxData<double,NUMCOMPS,HOST> Vector;

PROTO_KERNEL_START
	unsigned int InitializeStateF(State& a_U,
	                              const double a_gamma)
	{

		double gamma = a_gamma;
		double rho = 0.0;
		double p = 0.0;
		double u = 0.0;
		double v = 0.0;
		double w = 0.0;
		double Bx = 0.0;
		double By = 0.0;
		double Bz = 0.0;

		cout << "Here" << endl;
		
		// //////Modifying parameters for constant solution/////
		rho = 1.0;
		p = 1.0;

		double e = p/(gamma-1.0) + rho*(u*u+v*v+w*w)/2.0 + (Bx*Bx+By*By+Bz*Bz)/8.0/PI;

		a_U(0) = rho; //rho
		a_U(1) = rho*u; //Momentum-x
		a_U(2) = rho*v; //Momentum-y
		a_U(3) = rho*w; //Momentum-z
		a_U(4) = e; //Energy
		a_U(5) = Bx; //Bx
		a_U(6) = By; //By
		a_U(7) = Bz; //Bz


		return 0;
	}
	PROTO_KERNEL_END(InitializeStateF, InitializeState)


	void initializeState(BoxData<double,NUMCOMPS>& a_state,
	                     const double a_dx,
	                     const double a_dy,
	                     const double a_dz,
	                     const double a_gamma)
	{

		Box dbx0=a_state.box();
		Box dbx = dbx0.grow(NGHOST);
		Box dbx1 = dbx.grow(1);
		BoxData<double,NUMCOMPS,HOST> UBig(dbx1);
		forallInPlace(InitializeState,dbx1,UBig,a_gamma);
		// UBig = forall<double,NUMCOMPS>(InitializeState,a_gamma);
		Stencil<double> Lap2nd = Stencil<double>::Laplacian();
		Vector Lap = Lap2nd(UBig,dbx,1.0/24.0);
		UBig +=  Lap;
		a_state.setVal(0.0);
		a_state+=UBig;
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
			const char* varnames[8];
			varnames[0] = "density";
			varnames[1] = "Vx";
			varnames[2] = "Vy";
			varnames[3] = "Vz";
			varnames[4] = "p";
			varnames[5] = "Bx";
			varnames[6] = "By";
			varnames[7] = "Bz";

			double origin[DIM];
			for (int ii = 0; ii < DIM; ii++)
			{
				origin[ii] = 0.0;
			}
			WriteBoxData(filename_data.c_str(),out_data[*dit],varnames,origin,dx);
		}


	}


void WriteBoxData_array_nocoord(BoxData<double,NUMCOMPS>& out_data,
	                                       const double dx,
	                                       const double dy,
	                                       const double dz,
	                                       const string& filename_data)
	{


			const char* varnames[8];
			varnames[0] = "density";
			varnames[1] = "Vx";
			varnames[2] = "Vy";
			varnames[3] = "Vz";
			varnames[4] = "p";
			varnames[5] = "Bx";
			varnames[6] = "By";
			varnames[7] = "Bz";


			double origin[DIM];
			for (int ii = 0; ii < DIM; ii++)
			{
				origin[ii] = 0.0;
			}
			WriteBoxData(filename_data.c_str(),out_data,varnames,origin,dx);
	}	

int main(int argc, char* argv[])
{


#ifdef PR_MPI
	MPI_Init(&argc,&argv);
#endif

	//have to do this to get a time table
	PR_TIMER_SETFILE("proto.time.table");

	//PR_TIME("main");

	int pid = procID();

	// Defining inputs here until parmparse (or something similar becomes available)
	grid_type_global = 0;  // 0: 2D-Rectangular/3D-Rectangular;  1: 2D-Wavy/3D-Not implemented;  2: 2D-Polar/3D-Spherical
	pole_correction = false;
	double tstop = 0.5, CFL  = 0.01, domsizex = 1.0, domsizey = 1.0, domsizez = 1.0, gamma = 5.0/3.0;
	domainSizex = 64, domainSizey = 64, domainSizez = 64;
	BoxSize = 64;

	Box domain(Point::Zeros(),Point(domainSizex-1, domainSizey-1, domainSizez-1));

	array<bool,DIM> per;
	for(int idir = 0; idir < DIM; idir++){
		per[idir]=true;
	}

	ProblemDomain pd(domain,per);
	double dx = domsizex/domainSizex, dy = domsizey/domainSizey, dz = domsizez/domainSizez;

	MHDLevelDataState state(pd,BoxSize*Point::Ones(),dx, dy, dz, gamma);
	// InitializeMHDLevelDataState(state);
	(state.m_U).setToZero();
	
	for(DataIterator dit=state.m_U.begin(); *dit!=dit.end(); ++dit){
	initializeState((state.m_U)[*dit] ,state.m_dx, state.m_dy, state.m_dz,state.m_gamma);
	auto& state_i = state.m_U[*dit];
	WriteBoxData_array_nocoord(state_i, state.m_dx, state.m_dy, state.m_dz, "m_U");
	}
	
	LevelBoxData<double,NUMCOMPS> OUT[1];
	OUT[0].define(DisjointBoxLayout(pd,Point(domainSizex, domainSizey,  domainSizez)), Point::Zeros());
	(state.m_U).copyTo(OUT[0]);
	OUT[0].exchange();
	std::string filename_data="Output_"+std::to_string(0);
	WriteSinglePatchLevelData_nocoord(OUT[0], dx,dy,dz, filename_data);
	if(pid==0) cout << "Written .vtk file after step "<< 0 << endl;

	PR_TIMER_REPORT();

#ifdef PR_MPI
	MPI_Finalize();
#endif
}
