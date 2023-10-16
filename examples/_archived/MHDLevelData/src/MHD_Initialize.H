#pragma once

#include "Proto.H"
#include "MHDOp.H"

using namespace Proto;

namespace MHD_Initialize
{
	void initializeState(BoxData<double,NUMCOMPS>& a_state,
	                     const double a_dx,
	                     const double a_dy,
	                     const double a_dz,
	                     const double a_gamma);

	void initializeState_Spherical(BoxData<double,NUMCOMPS>& a_state,
						 BoxData<double,DIM*DIM>& a_detAA_avg,
						 BoxData<double,DIM*DIM>& a_detAA_inv_avg,
	                     BoxData<double,1>& a_r2rdot_avg,
	                     BoxData<double,1>& a_detA_avg,
	                     const double a_dx,
	                     const double a_dy,
	                     const double a_dz,
	                     const double a_gamma);		

	void InitializeStatecalc(BoxData<double,NUMCOMPS>& a_U,
	                         const BoxData<double,DIM>& a_x,
	                         const BoxData<double,DIM>& a_eta,
	                         const double a_gamma, 
							 const Box& a_dbx);
							 
	// PROTO_KERNEL_START
	// unsigned int InitializeStateF(const Point& a_pt,
	// 							  State& a_U,
	//                               const V& a_x,
	//                               const V& a_eta,
	//                               const double a_gamma);					 		 
}
