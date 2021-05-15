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
//#include "ParmParse.H"

extern double time_globalll;
extern bool limiter_apply;
extern bool slope_flattening_apply;
extern int Riemann_solver_type;

typedef BoxData<double,1,1,1> Scalar;
typedef BoxData<double,NUMCOMPS,1,1> Vector;

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

	#if DIM == 1
		rho = a_W(0);
		u   = a_W(1);
		p   = a_W(2);
		Bx  = a_W(3);
#endif
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
	                        const State& a_dot_pro_sum)
	{
		for (int i=0; i< NUMCOMPS; i++) {
			a_F_f_mapped1D(i) = a_N_s_d_ave_f(0)*a_F_ave_f(i) + a_dot_pro_sum(i)/12.0;
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
	          const BoxData<double,NUMCOMPS>& a_U,
	          const Box& a_rangeBox,
	          const double a_dx,
	          const double a_gamma,
	          Reduction<double>& a_Rxn,
	          BoxData<double,1>& Jacobian_ave,
	          bool a_computeMaxWaveSpeed,
	          bool a_callBCs)
	{
		// double nxq;
		//ParmParse pp;
		//pp.get("nx"                 , nxq);
		//cout << "nxq = " << nxq << endl;

		//cout << time_globalll << endl;
		Box dbx0 = a_U.box();
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
				m_interp_H[dir] = Stencil<double>::CellToEdgeH(dir);
				m_interp_L[dir] = Stencil<double>::CellToEdgeL(dir);
				m_interp_edge[dir] = Stencil<double>::CellToEdge(dir);
				m_divergence[dir] = Stencil<double>::FluxDivergence(dir);
				m_derivative[dir] = Stencil<double>::Derivative(1,dir,2);
			}
			initialized =  true;
		}

		using namespace std;
		a_Rhs.setVal(0.0);
		double gamma = a_gamma;

		if(a_callBCs)
		{
			BoxData<double, NUMCOMPS>& castU = const_cast<BoxData<double, NUMCOMPS> &>(a_U);
			int nghost = a_rangeBox.low()[0] - castU.box().low()[0];
			for(int idir = 0; idir < DIM; idir++)
			{
				protocommon::enforceSGBoundaryConditions(castU, nghost, idir);
			}
		}


		Vector a_U_demapped(dbx0);
		MHD_Mapping::JU_to_U_calc(a_U_demapped, a_U, Jacobian_ave, dbx0);


		Vector W_bar = forall<double,NUMCOMPS>(consToPrim,a_U_demapped, gamma);

		Vector U = m_deconvolve(a_U_demapped);

		Vector W    = forall<double,NUMCOMPS>(consToPrim,U, gamma);



		if(a_computeMaxWaveSpeed)
		{
			Scalar umax = forall<double>(waveSpeedBound,a_rangeBox,W, gamma);

			umax.absMax(a_Rxn); // returns 0 when used with CUDA
		}


		Vector W_ave = m_laplacian(W_bar,1.0/24.0);
		W_ave += W;

		// std::string filename="W_ave";
		// BoxData<double> W_ave_out = slice(W_ave,3);
		// WriteBoxData(filename.c_str(),W_ave_out,a_dx);

		for (int d = 0; d < DIM; d++)
		{

			//Vector W_ave_low = m_interp_L[d](W_ave);
			//Vector W_ave_high = m_interp_H[d](W_ave);

			Vector W_ave_low = m_interp_edge[d](W_ave);
			Vector W_ave_high = m_interp_edge[d](W_ave);

			Vector W_ave_low_lim_flat(dbx0),W_ave_high_lim_flat(dbx0);
			if (limiter_apply || slope_flattening_apply) {

				MHD_Limiters::MHD_Limiters(W_ave_low_lim_flat,W_ave_high_lim_flat,W_ave_low,W_ave_high,W_ave,W_bar,d,a_dx);

				W_ave_low = m_copy(W_ave_low_lim_flat);
				W_ave_high = m_copy(W_ave_high_lim_flat);
			}

#if DIM>1
			Vector W_low = m_deconvolve_f[d](W_ave_low);
			Vector W_high = m_deconvolve_f[d](W_ave_high);
#else
			Vector W_low = W_ave_low;
			Vector W_high = W_ave_high;
#endif

			Vector F_low(dbx0), F_high(dbx0), F_f(dbx0), F_ave_f(dbx0);
			Scalar Lambda_f(dbx0);
			Scalar N_s_d_ave_f(dbx0);
			Vector F_f_mapped(dbx0);
			F_f_mapped.setVal(0.0);
			for (int s = 0; s < DIM; s++) {
				if (Riemann_solver_type == 1) {
					MHD_Riemann_Solvers::Rusanov_Solver(F_f,W_low,W_high,s,gamma);
				}
				if (Riemann_solver_type == 2) {
					MHD_Riemann_Solvers::Roe8Wave_Solver(F_f,W_low,W_high,s,gamma);
				}
				MHD_Mapping::N_ave_f_calc_func(N_s_d_ave_f,s,d,a_dx);
#if DIM>1
				F_ave_f = m_convolve_f[d](F_f);
#else
				F_ave_f = F_f;
#endif
				Vector dot_pro_sum(dbx0);
				dot_pro_sum.setVal(0.0);
				for (int s_temp = 0; s_temp < DIM; s_temp++) {
					if (s_temp != d) {
						Scalar d_perp_N_s = m_derivative[s_temp](N_s_d_ave_f);
						Vector d_perp_F = m_derivative[s_temp](F_ave_f);
						Vector dot_pro = forall<double,NUMCOMPS>(dot_pro_calcF,d_perp_N_s,d_perp_F);
						dot_pro_sum += dot_pro;
					}
				}

				Vector F_f_mapped1D = forall<double,NUMCOMPS>(F_f_mapped1D_calc,F_ave_f,N_s_d_ave_f,dot_pro_sum);

				F_f_mapped += F_f_mapped1D;
			}
			a_Rhs += m_divergence[d](F_f_mapped);


			// if (d==0)
			// {
			// filename="RHS";
			// BoxData<double> RHS_out = slice(a_Rhs,2);
			// WriteBoxData(filename.c_str(),RHS_out,a_dx);
			// }

		}
		PR_TIME("EulerOp::operator::RHS*=-1.0/dx");
		a_Rhs *= -1./a_dx;
	}

}
