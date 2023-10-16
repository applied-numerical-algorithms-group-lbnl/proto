#include "Proto.H"
#include "MHD_Mapping.H"
#include "MHDOp.H"
#include "MHD_Initialize.H"
#include "MHD_Output_Writer.H"
#include "MHD_Input_Parsing.H"
extern Parsefrominputs inputs;

typedef BoxData<double,1,HOST> Scalar;
typedef BoxData<double,NUMCOMPS,HOST> Vector;

namespace MHD_Initialize {


	PROTO_KERNEL_START
	void dot_pro_calcFF(State& a_dot_pro,
	                    const Var<double,1>& a_d_perp_N_s,
	                    const State& a_d_perp_F)
	{
		for (int i=0; i< NUMCOMPS; i++) {
			a_dot_pro(i) = (a_d_perp_N_s(0)*a_d_perp_F(i));
			// a_dot_pro(i) = (a_d_perp_F(i));
		}
	}
	PROTO_KERNEL_END(dot_pro_calcFF, dot_pro_calcF)


	PROTO_KERNEL_START
	void InitializeStateF(const Point& a_pt,
								  State& a_U,
	                              V& a_x,
	                              V& a_eta,
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
		double x = a_x(0);
		double y = a_x(1);
#if DIM == 2
		double z = 0.0;
#endif
#if DIM == 3
		double z = a_x(2);
#endif
		// cout << a_pt [0] << " " << a_pt [1] << " " << a_pt [2] << endl;
		if (inputs.init_condition_type == 0) {
			// //////Modifying parameters for constant solution/////
			rho = 1.0;
			p = 1.0;
			
		}

		if (inputs.init_condition_type == 1) {
			// //////Modifying parameters for 2D current sheet problem/////
			rho = 1.0;
			p = 0.1;
			u = 0.1 * sin(2*PI*y);
			if (x >= 0.0 && x < 0.5) {By = 1.0;}
			if (x >= 0.5 && x < 1.5) {By = -1.0;}
			if (x >= 1.5 && x <= 2.0) {By = 1.0;}
		}

		if (inputs.init_condition_type == 2) {
			//////Modifying parameters for flow from 1 side///////
			rho = 1.0;
			p = 1.0;
			u = 1.0;
			v =  1.0;
		}

		if (inputs.init_condition_type == 3) {
			//////Modifying parameters for 2D Orszag Tang problem///////
			// Case 1:
			rho = gamma*((2.0 * gamma)/8.0/PI)*1.0;
			p = (2.0 * gamma)/8.0/M_PI;
			u = -sin(2.0*PI*y);
			v =  sin(2.0*PI*x);
			Bx = -sin(2.0*PI*y);
			By =  sin(4.0*PI*x);
		}

		if (inputs.init_condition_type == 4) {
			//////Modifying parameters for Alfven wave problem///////
			rho = 1.0;
			p = 1.0;
			u =  sin(2.0*PI*x);
			Bx = sin(2.0*PI*x);
		}

		if (inputs.init_condition_type == 5) {
			//////Modifying parameters for Acoustic pulse problem///////

			double rho_0 = 1.4;
			double delta_rho_0 = 0.14;
			double rad = sqrt((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5));
			if (rad < 0.5) {
				rho = rho_0 + delta_rho_0*exp(-16.0*rad*rad)*pow(cos(PI*rad),6.0);
			} else {
				rho = rho_0;
			}
			p = pow(rho/rho_0,gamma);
		}

		if (inputs.init_condition_type == 6) {
			//////Modifying parameters for Acoustic pulse problem in polar grid///////
			double pulsecenter_x = 0.0;
			double pulsecenter_y = 0.0;
			double rho_0 = 1.4;
			double delta_rho_0 = 0.04;
			double rad = sqrt((x-pulsecenter_x)*(x-pulsecenter_x) + (y-pulsecenter_y)*(y-pulsecenter_y));
			//if (rad < 0.6 && rad > 0.4){
			rho = rho_0 + delta_rho_0*exp(-400.0*(rad-0.5)*(rad-0.5))*pow(cos(PI*(rad-0.5)),16.0);
			//} else {
			//rho = rho_0;
			//}
			p = pow(rho/rho_0,gamma);
		}

		if (inputs.init_condition_type == 7) {
			//////Modifying parameters for cylindrical pulse problem in polar grid///////
			double rho_0 = 1.4;
			double delta_rho_0 = 0.14;
			if (x < 0.6 && x > 0.4) {
				rho = rho_0 + delta_rho_0*exp(-400.0*(x-0.5)*(x-0.5))*pow(cos(PI*(x-0.5)),16.0);
			} else {
				rho = rho_0;
			}
			p = pow(rho/rho_0,gamma);
		}

		if (inputs.init_condition_type == 8) {
			//////Modifying parameters for Shifted Acoustic pulse problem///////
			double pulsecenter_x = 0.4;
			double pulsecenter_y = 0.4;
			double rho_0 = 1.4;
			double delta_rho_0 = 0.14;
			double rad = sqrt((x-pulsecenter_x)*(x-pulsecenter_x) + (y-pulsecenter_y)*(y-pulsecenter_y));
			if (rad < 0.1) {
				rho = rho_0 + delta_rho_0*exp(-400.0*(rad)*(rad))*pow(cos(PI*(rad)),16.0);
			} else {
				rho = rho_0;
			}
			p = pow(rho/rho_0,gamma);
		}

		if (inputs.init_condition_type == 9) {
			//////Modifying parameters for Euler problem///////
			double rho0 = gamma;
			double p0 = 1.;
			rho = rho0;
			rho += .01*rho0*sin(2.0*2.0*PI*x);
			p = p0*pow(rho/rho0,gamma);
			double c0 = sqrt(gamma*p0/rho0);
			double c = sqrt(gamma*p/rho);
			u = 2*(c-c0)/(gamma-1.);
		}

		if (inputs.init_condition_type == 10) {
			//////Modifying parameters for radially out flow///////
			double theta = atan2(y,x);
			double rad = sqrt(y*y+x*x);
			rho = 1.0/rad/rad;
			p = rho;
			u = 1.0*cos(theta);
			v = 1.0*sin(theta);
			Bx = 1.0*cos(theta)/rad/rad;
			By = 1.0*sin(theta)/rad/rad;
		}

		if (inputs.init_condition_type == 11) {
			//////Modifying parameters for cylindrical pulse problem in cartesian grid///////
			double rho_0 = 1.4;
			double delta_rho_0 = 0.14;
			rho = rho_0 + delta_rho_0*exp(-400.0*(x-0.5)*(x-0.5))*pow(cos(PI*(x-0.5)),16.0);
			p = pow(rho/rho_0,gamma);
		}

		if (inputs.init_condition_type == 12) {
			//////Modifying parameters for Acoustic pulse problem///////

			double rho_0 = 1.4;
			double delta_rho_0 = 0.14;
			double rad = sqrt((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) + (z-0.5)*(z-0.5));
			if (rad < 0.5) {
				rho = rho_0 + delta_rho_0*exp(-16.0*rad*rad)*pow(cos(PI*rad),6.0);
			} else {
				rho = rho_0;
			}
			p = pow(rho/rho_0,gamma);
		}

		if (inputs.init_condition_type == 13) {
			//////Modifying parameters for 3D MHD blast wave///////
			rho = 1.0;
			double rad = sqrt((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) + (z-0.5)*(z-0.5));
			if (rad < 0.1) {
				p = 20.0;
			} else {
				p = 1.0;
			}
			Bx = 100.0/sqrt(4.0*PI);
		}

		if (inputs.init_condition_type == 14) {
			//////Modifying parameters for 2D MHD blast wave///////
			rho = 1.0;
			double rad = sqrt((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5));
			if (rad < 0.1) {
				p = 20.0;
			} else {
				p = 1.0;
			}
			Bx = 100.0/sqrt(4.0*PI);
		}

		if (inputs.init_condition_type == 15) {
			//////Modifying parameters for Acoustic pulse problem with Bx///////

			double rho_0 = 1.4;
			double delta_rho_0 = 0.14;
			double rad = sqrt((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5));
			if (rad < 0.5) {
				rho = rho_0 + delta_rho_0*exp(-16.0*rad*rad)*pow(cos(PI*rad),6.0);
			} else {
				rho = rho_0;
			}
			p = pow(rho/rho_0,gamma);
			Bx = 10.0/sqrt(4.0*PI);
		}

		if (inputs.init_condition_type == 16) {
			//////Modifying parameters for Acoustic pulse problem in Spherical grid///////
			double pulsecenter_x = 0.0;
			double pulsecenter_y = 0.0;
			double pulsecenter_z = 0.0;
			double rho_0 = 1.4;
			double delta_rho_0 = 0.14;
			double rad = sqrt((x-pulsecenter_x)*(x-pulsecenter_x) + (y-pulsecenter_y)*(y-pulsecenter_y) + (z-pulsecenter_z)*(z-pulsecenter_z));
			//if (rad < 0.6 && rad > 0.4){
			rho = rho_0 + delta_rho_0*exp(-400.0*(rad-0.5)*(rad-0.5))*pow(cos(PI*(rad-0.5)),16.0);
			//} else {
			//rho = rho_0;
			//}
			p = pow(rho/rho_0,gamma);
		}

		if (inputs.init_condition_type == 17) {
			//////Modifying parameters for Shifted Acoustic pulse problem in Spherical grid///////
			double rho_0 = 1.4;
			double delta_rho_0 = 0.14;

			double pulsecenter_x = 0.05;
			double pulsecenter_y = 0.0;
			double pulsecenter_z = 0.4;
			double rad = sqrt((x-pulsecenter_x)*(x-pulsecenter_x) + (y-pulsecenter_y)*(y-pulsecenter_y) + (z-pulsecenter_z)*(z-pulsecenter_z));

			rho = rho_0;

			if (rad < 0.1) {
				rho = rho_0 + delta_rho_0*exp(-400.0*(rad)*(rad))*pow(cos(PI*(rad)),16.0);
			}

			pulsecenter_x = 0.4;
			pulsecenter_y = 0.0;
			pulsecenter_z = 0.0;
			rad = sqrt((x-pulsecenter_x)*(x-pulsecenter_x) + (y-pulsecenter_y)*(y-pulsecenter_y) + (z-pulsecenter_z)*(z-pulsecenter_z));

			if (rad < 0.1) {
				rho = rho_0 + delta_rho_0*exp(-400.0*(rad)*(rad))*pow(cos(PI*(rad)),16.0);
			}

			p = pow(rho/rho_0,gamma);
		}

		if (inputs.init_condition_type == 18) {
			//////Modifying parameters for a spherical test///////
			double phi = atan2(y,x);
			rho = sin(phi);
			p = rho;
		}

		if (inputs.init_condition_type == 19) {
			//////Modifying parameters for velocity pulse problem in polar grid///////
			double rho_0 = 1.4;
			double delta_rho_0 = 0.4;
			double rad = sqrt((x)*(x) + (y)*(y));
			double theta = atan2(y,x);
			rho = rho_0;
			p = rho;
			//u = delta_rho_0*exp(-400.0*(rad-0.5)*(rad-0.5))*pow(cos(PI*(rad-0.5)),16.0)*cos(theta);
			//v = delta_rho_0*exp(-400.0*(rad-0.5)*(rad-0.5))*pow(cos(PI*(rad-0.5)),16.0)*sin(theta);

			if (abs(rad-0.5) < 0.1){
				double ang = 10.0*PI*(rad-0.4) - PI/2.0;
				u = (1.0*delta_rho_0+delta_rho_0*sin(ang))*cos(theta);
				v = (1.0*delta_rho_0+delta_rho_0*sin(ang))*sin(theta);
			} else {
				u = (0.0*delta_rho_0)*cos(theta);
				v = (0.0*delta_rho_0)*sin(theta);
			}
		}


		if (inputs.init_condition_type == 20) {
			//////Modifying parameters for theta dependent radially out flow ///////
			double theta = atan2(y,x);
			double rad = sqrt(y*y+x*x);
			rho = 1.0/rad/rad;
			p = rho;
			u = (5.0+abs(sin(theta)))*cos(theta);
			v = (5.0+abs(sin(theta)))*sin(theta);
			Bx = (2.0+abs(sin(theta)))*cos(theta)/rad/rad;
			By = (2.0+abs(sin(theta)))*sin(theta)/rad/rad;
		}

		if (inputs.init_condition_type == 21) {
			//////Modifying parameters for radially out flow in spherical grid///////

			double rad = sqrt(x*x+y*y+z*z);
			double phi = atan2(y,x);
			double theta = acos(z/rad);
			rho = 10.0*inputs.r_in*inputs.r_in/rad/rad;
			// rho = 10.0;
			p = 10.0*pow(inputs.r_in/rad,2.0*a_gamma);
			// p = 10.0;
			if (inputs.initialize_in_spherical_coords == 1){
				u = 5.0;
				v = 0.1*sin(theta)*cos(theta);
				w = 0.1;
			} else {
				u = 5.0*sin(theta)*cos(phi);
				v = 5.0*sin(theta)*sin(phi);
				w = 5.0*cos(theta);
			}
		}

		double e = p/(gamma-1.0) + rho*(u*u+v*v+w*w)/2.0 + (Bx*Bx+By*By+Bz*Bz)/8.0/PI;


#if DIM == 2
		a_U(0) = rho; //rho
		a_U(1) = rho*u; //Momentum-x
		a_U(2) = rho*v; //Momentum-y
		a_U(3) = e; //Energy
		a_U(4) = Bx; //Bx
		a_U(5) = By; //By
#endif
#if DIM == 3
		a_U(0) = rho; //rho
		a_U(1) = rho*u; //Momentum-x
		a_U(2) = rho*v; //Momentum-y
		a_U(3) = rho*w; //Momentum-z
		a_U(4) = e; //Energy
		a_U(5) = Bx; //Bx
		a_U(6) = By; //By
		a_U(7) = Bz; //Bz
#endif

		// return 0;
	}
	PROTO_KERNEL_END(InitializeStateF, InitializeState)


	void InitializeStatecalc(BoxData<double,NUMCOMPS>& a_U,
	                         const BoxData<double,DIM>& a_x,
	                         const BoxData<double,DIM>& a_eta,
	                         const double a_gamma, 
							 const Box& a_dbx)
	{	
		double gamma = a_gamma;
		// a_U = forall<double,NUMCOMPS>(InitializeState,a_x,a_eta, gamma);
		forallInPlace_p(InitializeState,a_U,a_x,a_eta,a_gamma);
	}



	void initializeState(BoxData<double,NUMCOMPS>& a_state,
	                     const double a_dx,
	                     const double a_dy,
	                     const double a_dz,
	                     const double a_gamma)
	{

		static Stencil<double> m_derivative[DIM];
		static Stencil<double> m_divergence[DIM];
		static Stencil<double> m_ahead_shift[DIM];
		static Stencil<double> m_behind_shift[DIM];
		static bool initialized = false;
		if(!initialized) {
			for (int dir = 0; dir < DIM; dir++) {
				m_derivative[dir] = Stencil<double>::Derivative(1,dir,2);
				m_divergence[dir] = Stencil<double>::FluxDivergence(dir);
				m_ahead_shift[dir] = 1.0*Shift(Point::Basis(dir)*(1));
				m_behind_shift[dir] = 1.0*Shift(Point::Basis(dir)*(-1));
			}
			initialized =  true;
		}

		Box dbx0=a_state.box();
		Box dbx = dbx0.grow(NGHOST);
		Box dbx1 = dbx.grow(1);
		BoxData<double,NUMCOMPS> UBig(dbx1);
		BoxData<double,DIM> eta(dbx1);
		MHD_Mapping::eta_calc(eta,dbx1,a_dx, a_dy, a_dz);
		BoxData<double,DIM> x(dbx1);
		
		MHD_Mapping::eta_to_x_calc(x,eta);
		
		// forallInPlace(InitializeState,dbx1,UBig,x,eta,a_gamma);
		// cout << "Here1" << endl;
		forallInPlace_p(InitializeState,UBig,x,eta,a_gamma);
		// cout << "Here2" << endl;
		{

		
		Stencil<double> Lap2nd = Stencil<double>::Laplacian();
		Vector Lap = Lap2nd(UBig,dbx,1.0/24.0);
		UBig +=  Lap;
		}
		// MHD_Output_Writer::WriteBoxData_array_nocoord(UBig, a_dx, a_dy, a_dz, "UBig");

		Scalar Jacobian_ave(dbx);

		if (inputs.grid_type_global == 2){
			MHD_Mapping::Jacobian_ave_sph_calc_func(Jacobian_ave, a_dx, a_dy, a_dz);
		} else {
			MHD_Mapping::Jacobian_Ave_calc(Jacobian_ave,a_dx, a_dy, a_dz, dbx);
		}

		a_state = forall<double,NUMCOMPS>(dot_pro_calcF, Jacobian_ave, UBig);
		for (int d=0; d<DIM; d++) {
			Vector d_UBig = m_derivative[d](UBig);
			Scalar d_Jacobian_ave = m_derivative[d](Jacobian_ave);
			Vector dot_pro = forall<double,NUMCOMPS>(dot_pro_calcF,d_Jacobian_ave,d_UBig);
			dot_pro *= 1./12.;
			a_state += dot_pro;
		}
		// MHD_Output_Writer::WriteBoxData_array_nocoord(a_state, a_dx, a_dy, a_dz, "a_state");
		// MHD_Output_Writer::WriteBoxData_array_nocoord(a_state, a_dx, a_dy, a_dz, "a_state2");
		// cout << "Here_in" << endl;
	}


	void initializeState_Spherical(BoxData<double,NUMCOMPS>& a_state,
						 BoxData<double,DIM*DIM>& a_detAA_avg,
						 BoxData<double,DIM*DIM>& a_detAA_inv_avg,
	                     BoxData<double,1>& a_r2rdot_avg,
	                     BoxData<double,1>& a_detA_avg,
	                     const double a_dx,
	                     const double a_dy,
	                     const double a_dz,
	                     const double a_gamma)
	{

		static Stencil<double> m_derivative[DIM];
		static Stencil<double> m_divergence[DIM];
		static Stencil<double> m_ahead_shift[DIM];
		static Stencil<double> m_behind_shift[DIM];
		static bool initialized = false;
		if(!initialized) {
			for (int dir = 0; dir < DIM; dir++) {
				m_derivative[dir] = Stencil<double>::Derivative(1,dir,2);
				m_divergence[dir] = Stencil<double>::FluxDivergence(dir);
				m_ahead_shift[dir] = 1.0*Shift(Point::Basis(dir)*(1));
				m_behind_shift[dir] = 1.0*Shift(Point::Basis(dir)*(-1));
			}
			initialized =  true;
		}

		Box dbx0=a_state.box();
		Box dbx = dbx0.grow(NGHOST);
		Box dbx1 = dbx.grow(1);
		BoxData<double,NUMCOMPS> UBig(dbx1);
		BoxData<double,DIM> eta(dbx1);
		MHD_Mapping::eta_calc(eta,dbx1,a_dx, a_dy, a_dz);
		BoxData<double,DIM> x(dbx1);		
		MHD_Mapping::eta_to_x_calc(x,eta);
		forallInPlace_p(InitializeState,UBig,x,eta,a_gamma);
		Stencil<double> Lap2nd = Stencil<double>::Laplacian();
		Vector Lap = Lap2nd(UBig,dbx,1.0/24.0);
		UBig +=  Lap;

		// MHD_Output_Writer::WriteBoxData_array_nocoord(UBig, a_dx, a_dy, a_dz, "UBig");
		MHD_Mapping::U_Sph_ave_to_JU_calc_func(a_state, UBig, a_detAA_avg, a_r2rdot_avg, a_detA_avg);
		// MHD_Output_Writer::WriteBoxData_array_nocoord(a_state, a_dx, a_dy, a_dz, "a_state");
		// MHD_Mapping::JU_to_U_Sph_ave_calc_func(UBig, a_state, a_detAA_inv_avg, a_r2rdot_avg, a_detA_avg, false);
		// MHD_Output_Writer::WriteBoxData_array_nocoord(UBig, a_dx, a_dy, a_dz, "UBig_again");

	}

}
