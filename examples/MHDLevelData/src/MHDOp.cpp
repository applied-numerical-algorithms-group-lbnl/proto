#include "Proto.H"
#include "MHDOp.H"
#include "CommonTemplates.H"
#include "Proto_Timer.H"
#include "Proto_WriteBoxData.H"
// For Chrono Timer (Talwinder)
#include <chrono>
#include <iostream>
#include <iomanip>

//////////////////////////////
#include "MHD_Limiters.H"
#include "MHD_Mapping.H"
#include "MHD_Riemann_Solvers.H"
#include "MHD_Output_Writer.H"
#include "MHD_Input_Parsing.H"

extern Parsefrominputs inputs;

typedef BoxData<double,1,HOST> Scalar;
typedef BoxData<double,NUMCOMPS,HOST> Vector;

namespace MHDOp {

	PROTO_KERNEL_START
	void
	consToPrimF(State&         a_W,
	            const State&   a_U,
	            double a_gamma)
	{
		double rho = a_U(0);
		double v2 = 0.0;
		double B2 = 0.0;
		double gamma = a_gamma;
		a_W(0) = rho;

		for (int i = 1; i <= DIM; i++)
		{
			double v, B;
			v = a_U(i) / rho;
			B = a_U(DIM+1+i);
			a_W(i) = v;
			a_W(DIM+1+i) = a_U(DIM+1+i);
			v2 += v*v;
			B2 += B*B;
		}

		a_W(NUMCOMPS-1-DIM) = (a_U(NUMCOMPS-1-DIM) - .5 * rho * v2  - B2/8.0/PI) * (gamma - 1.0);

	}
	PROTO_KERNEL_END(consToPrimF, consToPrim)







	PROTO_KERNEL_START
	void waveSpeedBoundF(Var<double,1>& a_speed,
	                     const State& a_W,
	                     double a_gamma)
	{
		double gamma = a_gamma;
		double rho=0., u=0., v=0., w=0., p=0., Bx=0., By=0., Bz=0., ce, af, B_mag, Bdir, udir;

#if DIM == 2
		rho = a_W(0);
		u   = a_W(1);
		v   = a_W(2);
		p   = a_W(3);
		Bx  = a_W(4);
		By  = a_W(5);
#endif
#if DIM == 3
		rho = a_W(0);
		u   = a_W(1);
		v   = a_W(2);
		w   = a_W(3);
		p   = a_W(4);
		Bx  = a_W(5);
		By  = a_W(6);
		Bz  = a_W(7);
#endif
		a_speed(0) = 0.0;
		for (int dir = 0; dir< DIM; dir++) {
			if (dir == 0) {
				Bdir = Bx;
				udir = u;
			};
			if (dir == 1) {
				Bdir = By;
				udir = v;
			};
			if (dir == 2) {
				Bdir = Bz;
				udir = w;
			};

			ce = sqrt(gamma*p/rho);
			B_mag = sqrt(Bx*Bx+By*By+Bz*Bz);
			af = 0.5*(sqrt((ce*ce)+( B_mag*B_mag/(4.0*PI*rho) )+( abs(Bdir)*ce/sqrt(PI*rho) ))+
			          sqrt((ce*ce)+( B_mag*B_mag/(4.0*PI*rho) )-( abs(Bdir)*ce/sqrt(PI*rho) ))) + abs(udir);
			if (af > a_speed(0)) {a_speed(0) = af;}
		}
	}
	PROTO_KERNEL_END(waveSpeedBoundF, waveSpeedBound)



	PROTO_KERNEL_START
	void dot_pro_calcFF(State& a_dot_pro,
	                    const Var<double,1>& a_d_perp_N_s,
	                    const State& a_d_perp_F)
	{
		for (int i=0; i< NUMCOMPS; i++) {
			a_dot_pro(i) = (a_d_perp_N_s(0)*a_d_perp_F(i));
		}
	}
	PROTO_KERNEL_END(dot_pro_calcFF, dot_pro_calcF)



	PROTO_KERNEL_START
	void F_f_mapped1D_calcF(State& a_F_f_mapped1D,
	                        const State& a_F_ave_f,
	                        const Var<double,1>& a_N_s_d_ave_f,
	                        const State& a_dot_pro_sum,
	                        const double a_dx_d)
	{
		for (int i=0; i< NUMCOMPS; i++) {
			// a_F_f_mapped1D(i) = -(a_N_s_d_ave_f(0)*a_F_ave_f(i) + a_dot_pro_sum(i)/12.0)/a_dx_d;
			a_F_f_mapped1D(i) = (a_N_s_d_ave_f(0)*a_F_ave_f(i) + a_dot_pro_sum(i)/12.0);
		}
	}
	PROTO_KERNEL_END(F_f_mapped1D_calcF, F_f_mapped1D_calc)


	void consToPrimcalc(BoxData<double,NUMCOMPS>& W_bar,
	                    const BoxData<double,NUMCOMPS>& a_U_demapped,
	                    const double gamma)
	{
		W_bar = forall<double,NUMCOMPS>(consToPrim,a_U_demapped, gamma);
	}



	void step(BoxData<double,NUMCOMPS>& a_Rhs,
	          const BoxData<double,NUMCOMPS>& a_JU_ave,
	          const Box& a_rangeBox,
	          const double a_dx,
	          const double a_dy,
	          const double a_dz,
	          const double a_gamma,
	          Reduction<double>& a_Rxn,
	          BoxData<double,1>& a_Jacobian_ave,
	          BoxData<double,DIM*DIM>& a_N_ave_f,
	          bool a_computeMaxWaveSpeed,
	          bool a_callBCs)
	{

		Box dbx0 = a_JU_ave.box();
		//Box dbx1 = dbx0.grow(1-NGHOST);
		Box dbx1 = dbx0;
		Box dbx2 = dbx0.grow(0-NGHOST);
		static Stencil<double> m_laplacian;
		static Stencil<double> m_deconvolve;
		static Stencil<double> m_copy;
		static Stencil<double> m_laplacian_f[DIM];
		static Stencil<double> m_deconvolve_f[DIM];
		static Stencil<double> m_convolve_f[DIM];
		static Stencil<double> m_interp_H[DIM];
		static Stencil<double> m_interp_L[DIM];
		static Stencil<double> m_interp_edge[DIM];
		static Stencil<double> m_divergence[DIM];
		static Stencil<double> m_derivative[DIM];
		static bool initialized = false;
		if(!initialized)
		{
			m_laplacian = Stencil<double>::Laplacian();
			m_deconvolve = (-1.0/24.0)*m_laplacian + (1.0)*Shift(Point::Zeros());
			m_copy = 1.0*Shift(Point::Zeros());
			for (int dir = 0; dir < DIM; dir++)
			{
				m_laplacian_f[dir] = Stencil<double>::LaplacianFace(dir);
				m_deconvolve_f[dir] = (-1.0/24.0)*m_laplacian_f[dir] + 1.0*Shift(Point::Zeros());
				m_convolve_f[dir] = (1.0/24.0)*m_laplacian_f[dir] + 1.0*Shift(Point::Zeros());
				m_interp_H[dir] = Stencil<double>::CellToFaceH(dir);
				m_interp_L[dir] = Stencil<double>::CellToFaceL(dir);
				m_interp_edge[dir] = Stencil<double>::CellToFace(dir);
				m_divergence[dir] = Stencil<double>::FluxDivergence(dir);
				m_derivative[dir] = Stencil<double>::Derivative(1,dir,2);
			}
			initialized =  true;
		}

		using namespace std;
		a_Rhs.setVal(0.0);
		double gamma = a_gamma;
		double dxd[3] = {a_dx, a_dy, a_dz};
		if(a_callBCs)
		{
			BoxData<double, NUMCOMPS>& castU = const_cast<BoxData<double, NUMCOMPS> &>(a_JU_ave);
			int nghost = a_rangeBox.low()[0] - castU.box().low()[0];
			for(int idir = 0; idir < DIM; idir++)
			{
				protocommon::enforceSGBoundaryConditions(castU, nghost, idir);
			}
		}


		Vector a_U_ave(dbx0);
		MHD_Mapping::JU_to_U_calc(a_U_ave, a_JU_ave, a_Jacobian_ave, dbx0);
		Vector W_bar = forall<double,NUMCOMPS>(consToPrim,a_U_ave, gamma);
		Vector U = m_deconvolve(a_U_ave);
		Vector W  = forall<double,NUMCOMPS>(consToPrim,U, gamma);
		Vector W_ave = m_laplacian(W_bar,1.0/24.0);
		W_ave += W;


		// Vector a_JU_temp_ave(dbx0);
		// a_JU_temp_ave.setVal(0.0);
		// a_JU_temp_ave += a_JU_ave;
		// MHD_Output_Writer::WriteBoxData_array_nocoord(a_JU_temp_ave, a_dx, a_dy, a_dz, "a_JU_ave_m1");
		// MHD_Output_Writer::WriteBoxData_array_nocoord(W_ave, a_dx, a_dy, a_dz, "W_ave_m1");

		for (int d = 0; d < DIM; d++)
		{
			Vector W_ave_low_temp(dbx0), W_ave_high_temp(dbx0);
			Vector W_ave_low(dbx0), W_ave_high(dbx0);
			// W_ave_low_temp = m_interp_edge[d](W_ave);
			// W_ave_high_temp = m_copy(W_ave_low_temp);
			W_ave_low_temp = m_interp_L[d](W_ave);
			W_ave_high_temp = m_interp_H[d](W_ave);
			MHD_Limiters::MHD_Limiters(W_ave_low,W_ave_high,W_ave_low_temp,W_ave_high_temp,W_ave,W_bar,d,a_dx, a_dy, a_dz);

			Vector W_low = m_deconvolve_f[d](W_ave_low);
			Vector W_high = m_deconvolve_f[d](W_ave_high);


			Vector F_f(dbx1), F_ave_f(dbx1);
			Vector F_f_mapped(dbx1);
			// Vector F_f_mapped_noghost(dbx2);
			F_f_mapped.setVal(0.0);
			// F_f_mapped_noghost.setVal(0.0);
			double dx_d = dxd[d];


			for (int s = 0; s < DIM; s++) {
				if (inputs.Riemann_solver_type == 1) {
					MHD_Riemann_Solvers::Rusanov_Solver(F_f,W_low,W_high,s,gamma);
				}
				if (inputs.Riemann_solver_type == 2) {
					MHD_Riemann_Solvers::Roe8Wave_Solver(F_f,W_low,W_high,s,gamma);
				}
				Scalar N_s_d_ave_f = slice(a_N_ave_f,d*DIM+s);

				F_ave_f = m_convolve_f[d](F_f);

				Vector dot_pro_sum(dbx1);
				dot_pro_sum.setVal(0.0);
				for (int s_temp = 0; s_temp < DIM; s_temp++) {
					if (s_temp != d) {
						Scalar d_perp_N_s = m_derivative[s_temp](N_s_d_ave_f);
						Vector d_perp_F = m_derivative[s_temp](F_ave_f);
						Vector dot_pro = forall<double,NUMCOMPS>(dot_pro_calcF,d_perp_N_s,d_perp_F);
						dot_pro_sum += dot_pro;
					}
				}

				Vector F_f_mapped1D = forall<double,NUMCOMPS>(F_f_mapped1D_calc,F_ave_f,N_s_d_ave_f,dot_pro_sum,dx_d);


				F_f_mapped += F_f_mapped1D;

			}
			Vector Rhs_d = m_divergence[d](F_f_mapped);
			//PR_TIME("EulerOp::operator::RHS*=-1.0/dx");
			Rhs_d *= -1./dx_d;
			a_Rhs += Rhs_d;
			// F_f_mapped_noghost += F_f_mapped;
			// if (d==0) MHD_Output_Writer::WriteBoxData_array_nocoord(F_f_mapped_noghost, a_dx, a_dy, a_dz, "F0_m1");
			// if (d==1) MHD_Output_Writer::WriteBoxData_array_nocoord(F_f_mapped_noghost, a_dx, a_dy, a_dz, "F1_m1");
			// if (d==2) MHD_Output_Writer::WriteBoxData_array_nocoord(F_f_mapped_noghost, a_dx, a_dy, a_dz, "F2_m1");

		}
	}







	void step_spherical(BoxData<double,NUMCOMPS>& a_Rhs,
	                    const BoxData<double,NUMCOMPS>& a_JU_ave,
	                    const Box& a_rangeBox,
	                    const double a_dx,
	                    const double a_dy,
	                    const double a_dz,
	                    const double a_gamma,
	                    Reduction<double>& a_Rxn,
	                    BoxData<double,1>& a_Jacobian_ave,
	                    BoxData<double,DIM*DIM>& a_N_ave_f,
						BoxData<double,DIM*DIM>& a_A_1_avg,
						BoxData<double,DIM*DIM>& a_A_2_avg,
						BoxData<double,DIM*DIM>& a_A_3_avg,
	                    BoxData<double,DIM*DIM>& a_detAA_avg,
	                    BoxData<double,DIM*DIM>& a_detAA_inv_avg,
	                    BoxData<double,1>& a_r2rdot_avg,
	                    BoxData<double,1>& a_detA_avg,
	                    BoxData<double,1>& a_r2detA_1_avg,
	                    BoxData<double,DIM*DIM>& a_r2detAA_1_avg,
	                    BoxData<double,DIM>& a_r2detAn_1_avg,
	                    BoxData<double,1>& a_rrdotdetA_2_avg,
	                    BoxData<double,DIM*DIM>& a_rrdotdetAA_2_avg,
	                    BoxData<double,DIM>& a_rrdotd3ncn_2_avg,
	                    BoxData<double,1>& a_rrdotdetA_3_avg,
	                    BoxData<double,DIM*DIM>& a_rrdotdetAA_3_avg,
	                    BoxData<double,DIM>& a_rrdotncd2n_3_avg,
	                    bool a_computeMaxWaveSpeed,
	                    bool a_callBCs)
	{

		Box dbx0 = a_JU_ave.box();
		Box dbx1 = dbx0.grow(NGHOST-NGHOST);
		//Box dbx1 = dbx0;
		Box dbx2 = dbx0.grow(1-NGHOST);
		static Stencil<double> m_laplacian;
		static Stencil<double> m_deconvolve;
		static Stencil<double> m_copy;
		static Stencil<double> m_laplacian_f[DIM];
		static Stencil<double> m_deconvolve_f[DIM];
		static Stencil<double> m_convolve_f[DIM];
		static Stencil<double> m_interp_H[DIM];
		static Stencil<double> m_interp_L[DIM];
		static Stencil<double> m_interp_edge[DIM];
		static Stencil<double> m_divergence[DIM];
		static Stencil<double> m_derivative[DIM];
		static bool initialized = false;
		if(!initialized)
		{
			m_laplacian = Stencil<double>::Laplacian();
			m_deconvolve = (-1.0/24.0)*m_laplacian + (1.0)*Shift(Point::Zeros());
			m_copy = 1.0*Shift(Point::Zeros());
			for (int dir = 0; dir < DIM; dir++)
			{
				m_laplacian_f[dir] = Stencil<double>::LaplacianFace(dir);
				m_deconvolve_f[dir] = (-1.0/24.0)*m_laplacian_f[dir] + 1.0*Shift(Point::Zeros());
				m_convolve_f[dir] = (1.0/24.0)*m_laplacian_f[dir] + 1.0*Shift(Point::Zeros());
				m_interp_H[dir] = Stencil<double>::CellToFaceH(dir);
				m_interp_L[dir] = Stencil<double>::CellToFaceL(dir);
				m_interp_edge[dir] = Stencil<double>::CellToFace(dir);
				m_divergence[dir] = Stencil<double>::FluxDivergence(dir);
				m_derivative[dir] = Stencil<double>::Derivative(1,dir,2);
			}
			initialized =  true;
		}

		using namespace std;
		a_Rhs.setVal(0.0);
		double gamma = a_gamma;
		double dxd[3] = {a_dx, a_dy, a_dz}; // Because now its r, theta, phi

		Vector a_U_Sph_ave(dbx0);
		MHD_Mapping::JU_to_U_Sph_ave_calc_func(a_U_Sph_ave, a_JU_ave, a_detAA_inv_avg, a_r2rdot_avg, a_detA_avg, false);
		Vector W_bar = forall<double,NUMCOMPS>(consToPrim,a_U_Sph_ave, gamma);
		Vector U = m_deconvolve(a_U_Sph_ave);
		Vector W  = forall<double,NUMCOMPS>(consToPrim,U, gamma);
		Vector W_ave = m_laplacian(W_bar,1.0/24.0);
		W_ave += W;

		for (int d = 0; d < DIM; d++)
		{
			Vector W_ave_low_temp(dbx0), W_ave_high_temp(dbx0);
			Vector W_ave_low(dbx0), W_ave_high(dbx0);
			Vector W_low_cart(dbx0), W_high_cart(dbx0);
			W_ave_low_temp = m_interp_L[d](W_ave);
			W_ave_high_temp = m_interp_H[d](W_ave);
			// W_ave_low_temp = m_interp_edge[d](W_ave);
			// W_ave_high_temp = m_copy(W_ave_low_temp);		
			MHD_Limiters::MHD_Limiters(W_ave_low,W_ave_high,W_ave_low_temp,W_ave_high_temp,W_ave,W_bar,d,a_dx, a_dy, a_dz);			
			Vector W_low = m_deconvolve_f[d](W_ave_low);
			Vector W_high = m_deconvolve_f[d](W_ave_high);
			Vector F_ave_f(dbx0), F_f(dbx0);
			F_ave_f.setVal(0.0);
			F_f.setVal(0.0);
			// Vector F_f_mapped_noghost(dbx2);
			// F_f_mapped_noghost.setVal(0.0);
			double dx_d = dxd[d];
			MHD_Riemann_Solvers::Spherical_Riemann_Solver(F_ave_f, W_low, W_high, W_low, W_high, W_ave_low, W_ave_high, a_r2detA_1_avg, a_r2detAA_1_avg, a_r2detAn_1_avg, a_rrdotdetA_2_avg, a_rrdotdetAA_2_avg, a_rrdotd3ncn_2_avg, a_rrdotdetA_3_avg, a_rrdotdetAA_3_avg, a_rrdotncd2n_3_avg, d, gamma, a_dx, a_dy, a_dz);	
			Vector Rhs_d = m_divergence[d](F_ave_f);
			//PR_TIME("EulerOp::operator::RHS*=-1.0/dx");
			Rhs_d *= -1./dx_d;
			a_Rhs += Rhs_d;
			// F_f_mapped_noghost += F_ave_f;
			// if (d==0) MHD_Output_Writer::WriteBoxData_array_nocoord(F_f_mapped_noghost, a_dx, a_dy, a_dz, "F0_m1");
			// if (d==1) MHD_Output_Writer::WriteBoxData_array_nocoord(F_f_mapped_noghost, a_dx, a_dy, a_dz, "F1_m1");
			// if (d==2) MHD_Output_Writer::WriteBoxData_array_nocoord(F_f_mapped_noghost, a_dx, a_dy, a_dz, "F2_m1");

		}
	}


}
