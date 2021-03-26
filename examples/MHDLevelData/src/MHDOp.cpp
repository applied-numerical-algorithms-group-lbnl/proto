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
//#include "ParmParse.H"

extern double time_globalll;

typedef BoxData<double,1,1,1> Scalar;
typedef BoxData<double,NUMCOMPS,1,1> Vector;

namespace MHDOp {
	
	
	
	PROTO_KERNEL_START
	void rusanovStateF(State& a_out,
						const State& a_W_lo,
						const State& a_W_hi,
						const State& a_F_lo,
						const State& a_F_hi,
						const Var<double,1>& a_lambda,
						int   a_dir,
						double a_gamma)
	{
		double gamma = a_gamma;
		
#if DIM == 1
		double rho_lo, rhou_lo, e_lo, p_lo, Bx_lo, B2_lo, v2_lo;
		double rho_hi, rhou_hi, e_hi, p_hi, Bx_hi, B2_hi, v2_hi;
		rho_lo  = a_W_lo(0);
		rho_hi  = a_W_hi(0);
		rhou_lo = rho_lo*a_W_lo(1);
		rhou_hi = rho_hi*a_W_hi(1);
		p_lo    = a_W_lo(2);
		p_hi    = a_W_hi(2);
		Bx_lo   = a_W_lo(3);
		Bx_hi   = a_W_hi(3);
		B2_lo   = Bx_lo*Bx_lo;
		B2_hi   = Bx_hi*Bx_hi;
		v2_lo   = a_W_lo(1)*a_W_lo(1);
		v2_hi   = a_W_hi(1)*a_W_hi(1);
		e_lo    = p_lo/(gamma-1.0) + rho_lo*v2_lo/2.0 + B2_lo/8.0/PI;
		e_hi    = p_hi/(gamma-1.0) + rho_hi*v2_hi/2.0 + B2_hi/8.0/PI;
		a_out(0) = 0.5*(a_F_hi(0) + a_F_lo(0) - abs(a_lambda(0))*(rho_hi-rho_lo));
		a_out(1) = 0.5*(a_F_hi(1) + a_F_lo(1) - abs(a_lambda(0))*(rhou_hi-rhou_lo));
		a_out(2) = 0.5*(a_F_hi(2) + a_F_lo(2) - abs(a_lambda(0))*(e_hi-e_lo));
		a_out(3) = 0.5*(a_F_hi(3) + a_F_lo(3) - abs(a_lambda(0))*(Bx_hi-Bx_lo));	
#endif    
		
#if DIM == 2
		double rho_lo, rhou_lo, rhov_lo, e_lo, p_lo, Bx_lo, By_lo, B2_lo, v2_lo;
		double rho_hi, rhou_hi, rhov_hi, e_hi, p_hi, Bx_hi, By_hi, B2_hi, v2_hi;
		rho_lo  = a_W_lo(0);
		rho_hi  = a_W_hi(0);
		rhou_lo = rho_lo*a_W_lo(1);
		rhou_hi = rho_hi*a_W_hi(1);
		rhov_lo = rho_lo*a_W_lo(2);
		rhov_hi = rho_hi*a_W_hi(2);
		p_lo    = a_W_lo(3);
		p_hi    = a_W_hi(3);
		Bx_lo   = a_W_lo(4);
		Bx_hi   = a_W_hi(4);
		By_lo   = a_W_lo(5);
		By_hi   = a_W_hi(5);
		B2_lo   = Bx_lo*Bx_lo + By_lo*By_lo;
		B2_hi   = Bx_hi*Bx_hi + By_hi*By_hi;
		v2_lo   = a_W_lo(1)*a_W_lo(1) + a_W_lo(2)*a_W_lo(2);
		v2_hi   = a_W_hi(1)*a_W_hi(1) + a_W_hi(2)*a_W_hi(2);
		e_lo    = p_lo/(gamma-1.0) + rho_lo*v2_lo/2.0 + B2_lo/8.0/PI;
		e_hi    = p_hi/(gamma-1.0) + rho_hi*v2_hi/2.0 + B2_hi/8.0/PI;
		a_out(0) = 0.5*(a_F_hi(0) + a_F_lo(0) - abs(a_lambda(0))*(rho_hi-rho_lo));
		a_out(1) = 0.5*(a_F_hi(1) + a_F_lo(1) - abs(a_lambda(0))*(rhou_hi-rhou_lo));
		a_out(2) = 0.5*(a_F_hi(2) + a_F_lo(2) - abs(a_lambda(0))*(rhov_hi-rhov_lo));
		a_out(3) = 0.5*(a_F_hi(3) + a_F_lo(3) - abs(a_lambda(0))*(e_hi-e_lo));
		a_out(4) = 0.5*(a_F_hi(4) + a_F_lo(4) - abs(a_lambda(0))*(Bx_hi-Bx_lo));
		a_out(5) = 0.5*(a_F_hi(5) + a_F_lo(5) - abs(a_lambda(0))*(By_hi-By_lo));
#endif 

#if DIM == 3
		double rho_lo, rhou_lo, rhov_lo, rhow_lo, e_lo, p_lo, Bx_lo, By_lo, Bz_lo, B2_lo, v2_lo;
		double rho_hi, rhou_hi, rhov_hi, rhow_hi, e_hi, p_hi, Bx_hi, By_hi, Bz_hi, B2_hi, v2_hi;
		rho_lo  = a_W_lo(0);
		rho_hi  = a_W_hi(0);
		rhou_lo = rho_lo*a_W_lo(1);
		rhou_hi = rho_hi*a_W_hi(1);
		rhov_lo = rho_lo*a_W_lo(2);
		rhov_hi = rho_hi*a_W_hi(2);
		rhow_lo = rho_lo*a_W_lo(3);
		rhow_hi = rho_hi*a_W_hi(3);
		p_lo    = a_W_lo(4);
		p_hi    = a_W_hi(4);
		Bx_lo   = a_W_lo(5);
		Bx_hi   = a_W_hi(5);
		By_lo   = a_W_lo(6);
		By_hi   = a_W_hi(6);
		Bz_lo   = a_W_lo(7);
		Bz_hi   = a_W_hi(7);
		B2_lo   = Bx_lo*Bx_lo + By_lo*By_lo + Bz_lo*Bz_lo;
		B2_hi   = Bx_hi*Bx_hi + By_hi*By_hi + Bz_hi*Bz_hi;
		v2_lo   = a_W_lo(1)*a_W_lo(1) + a_W_lo(2)*a_W_lo(2) + a_W_lo(3)*a_W_lo(3);
		v2_hi   = a_W_hi(1)*a_W_hi(1) + a_W_hi(2)*a_W_hi(2) + a_W_hi(3)*a_W_hi(3);
		e_lo    = p_lo/(gamma-1.0) + rho_lo*v2_lo/2.0 + B2_lo/8.0/PI;
		e_hi    = p_hi/(gamma-1.0) + rho_hi*v2_hi/2.0 + B2_hi/8.0/PI;
		a_out(0) = 0.5*(a_F_hi(0) + a_F_lo(0) - abs(a_lambda(0))*(rho_hi-rho_lo));
		a_out(1) = 0.5*(a_F_hi(1) + a_F_lo(1) - abs(a_lambda(0))*(rhou_hi-rhou_lo));
		a_out(2) = 0.5*(a_F_hi(2) + a_F_lo(2) - abs(a_lambda(0))*(rhov_hi-rhov_lo));
		a_out(3) = 0.5*(a_F_hi(3) + a_F_lo(3) - abs(a_lambda(0))*(rhow_hi-rhow_lo));
		a_out(4) = 0.5*(a_F_hi(4) + a_F_lo(4) - abs(a_lambda(0))*(e_hi-e_lo));
		a_out(5) = 0.5*(a_F_hi(5) + a_F_lo(5) - abs(a_lambda(0))*(Bx_hi-Bx_lo));
		a_out(6) = 0.5*(a_F_hi(6) + a_F_lo(6) - abs(a_lambda(0))*(By_hi-By_lo));
		a_out(7) = 0.5*(a_F_hi(7) + a_F_lo(7) - abs(a_lambda(0))*(Bz_hi-Bz_lo));
#endif 	
	}
	PROTO_KERNEL_END(rusanovStateF, rusanovState)
	
	
	PROTO_KERNEL_START
	void lambdacalcF(Var<double,1>& a_lambda,
					 const State& a_W_low,
					 const State& a_W_high,
					 int a_d,
					 double a_gamma)
	{
		double gamma = a_gamma;
		double rho=0., u=0., v=0., w=0., p=0., Bx=0., By=0., Bz=0., ce, af, B_mag, Bdir, udir;
#if DIM == 1
		rho = 0.5 * (a_W_low(0) + a_W_high(0));
		u   = 0.5 * (a_W_low(1) + a_W_high(1));
		p   = 0.5 * (a_W_low(2) + a_W_high(2));
		Bx  = 0.5 * (a_W_low(3) + a_W_high(3));
#endif	
#if DIM == 2
		rho = 0.5 * (a_W_low(0) + a_W_high(0));
		u   = 0.5 * (a_W_low(1) + a_W_high(1));
		v   = 0.5 * (a_W_low(2) + a_W_high(2));
		p   = 0.5 * (a_W_low(3) + a_W_high(3));
		Bx  = 0.5 * (a_W_low(4) + a_W_high(4));
		By  = 0.5 * (a_W_low(5) + a_W_high(5));
#endif	
#if DIM == 3
		rho = 0.5 * (a_W_low(0) + a_W_high(0));
		u   = 0.5 * (a_W_low(1) + a_W_high(1));
		v   = 0.5 * (a_W_low(2) + a_W_high(2));
		w   = 0.5 * (a_W_low(3) + a_W_high(3));
		p   = 0.5 * (a_W_low(4) + a_W_high(4));
		Bx  = 0.5 * (a_W_low(5) + a_W_high(5));
		By  = 0.5 * (a_W_low(6) + a_W_high(6));
		Bz  = 0.5 * (a_W_low(7) + a_W_high(7));
#endif	
		if (a_d == 0){
			Bdir = Bx;
			udir = u;
		};
		if (a_d == 1){
			Bdir = By;
			udir = v;
		};
		if (a_d == 2){
			Bdir = Bz;
			udir = w;
		};
		
		ce = sqrt(gamma*p/rho);
		B_mag = sqrt(Bx*Bx+By*By+Bz*Bz);
		af = 0.5*(sqrt((ce*ce)+( B_mag*B_mag/(4.0*PI*rho) )+( abs(Bdir)*ce/sqrt(PI*rho) ))+
			  sqrt((ce*ce)+( B_mag*B_mag/(4.0*PI*rho) )-( abs(Bdir)*ce/sqrt(PI*rho) )));
		a_lambda(0) = af + abs(udir);
	}
	PROTO_KERNEL_END(lambdacalcF, lambdacalc) 
	
	

  
  PROTO_KERNEL_START
	void 
	consToPrimF(State&         a_W, 
				const State&   a_U,
				double         a_gamma)
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
	void getFluxF(State& a_F, const State& a_W, 
				  int    a_d,
				  double a_gamma)
	  {
		double rho, p0, e, v2, B2, vB;  
		double gamma = a_gamma;
		double mult1=0.;
		double mult2=0.;
		double mult3=0.;	
		if (a_d == 0){mult1=1.0;}
		if (a_d == 1){mult2=1.0;}
		if (a_d == 2){mult3=1.0;}
		rho = a_W(0);
#if DIM==1	
		v2 = a_W(1)*a_W(1);
		B2 = a_W(3)*a_W(3); 
		vB = a_W(1)*a_W(3);
		p0 = a_W(2) + B2/8.0/PI;
		e  = a_W(2)/(gamma-1.0) + rho*v2/2.0 + B2/8.0/PI;
		a_F(0) = rho*a_W(1+a_d);
		a_F(1) = rho*a_W(1)*a_W(1+a_d) + mult1*p0 - a_W(3)*a_W(3+a_d)/4.0/PI;
		a_F(2) = (e+p0)*a_W(1+a_d) - a_W(3+a_d)*vB/4.0/PI;
		a_F(3) = 0.0;
#endif
		
#if DIM==2	
		v2 = a_W(1)*a_W(1) + a_W(2)*a_W(2);
		B2 = a_W(4)*a_W(4) + a_W(5)*a_W(5);
		vB = a_W(1)*a_W(4) + a_W(2)*a_W(5);
		p0 = a_W(3) + B2/8.0/PI;
		e  = a_W(3)/(gamma-1.0) + rho*v2/2.0 + B2/8.0/PI;
		a_F(0) = rho*a_W(1+a_d);
		a_F(1) = rho*a_W(1)*a_W(1+a_d) + mult1*p0 - a_W(4)*a_W(4+a_d)/4.0/PI;
		a_F(2) = rho*a_W(2)*a_W(1+a_d) + mult2*p0 - a_W(5)*a_W(4+a_d)/4.0/PI;
		a_F(3) = (e+p0)*a_W(1+a_d) - a_W(4+a_d)*vB/4.0/PI;
		a_F(4) = mult2*(a_W(1+a_d)*a_W(4) - a_W(1)*a_W(4+a_d));
		a_F(5) = mult1*(a_W(1+a_d)*a_W(5) - a_W(2)*a_W(4+a_d));	
#endif

#if DIM==3	
		v2 = a_W(1)*a_W(1) + a_W(2)*a_W(2) + a_W(3)*a_W(3);
		B2 = a_W(5)*a_W(5) + a_W(6)*a_W(6) + a_W(7)*a_W(7);
		vB = a_W(1)*a_W(5) + a_W(2)*a_W(6) + a_W(3)*a_W(7);
		p0 = a_W(4) + B2/8.0/PI;
		e  = a_W(4)/(gamma-1.0) + rho*v2/2.0 + B2/8.0/PI;
		a_F(0) = rho*a_W(1+a_d);
		a_F(1) = rho*a_W(1)*a_W(1+a_d) + mult1*p0 - a_W(5)*a_W(5+a_d)/4.0/PI;
		a_F(2) = rho*a_W(2)*a_W(1+a_d) + mult2*p0 - a_W(6)*a_W(5+a_d)/4.0/PI;
		a_F(3) = rho*a_W(3)*a_W(1+a_d) + mult3*p0 - a_W(7)*a_W(5+a_d)/4.0/PI;
		a_F(4) = (e+p0)*a_W(1+a_d) - a_W(5+a_d)*vB/4.0/PI;
		a_F(5) = (mult2+mult3)*(a_W(1+a_d)*a_W(5) - a_W(1)*a_W(5+a_d));
		a_F(6) = (mult1+mult3)*(a_W(1+a_d)*a_W(6) - a_W(2)*a_W(5+a_d));
		a_F(7) = (mult1+mult2)*(a_W(1+a_d)*a_W(7) - a_W(3)*a_W(5+a_d));
#endif	  
	  }
	  PROTO_KERNEL_END(getFluxF, getFlux)
  

    PROTO_KERNEL_START
	void waveSpeedBoundF(Var<double,1>& a_speed,
						 const State& a_W,
						 double       a_gamma)
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
		for (int dir = 0; dir< DIM; dir++){
			if (dir == 0){
				Bdir = Bx;
				udir = u;
			};
			if (dir == 1){
				Bdir = By;
				udir = v;
			};
			if (dir == 2){
				Bdir = Bz;
				udir = w;
			};
		
		ce = sqrt(gamma*p/rho);
		B_mag = sqrt(Bx*Bx+By*By+Bz*Bz);
		af = 0.5*(sqrt((ce*ce)+( B_mag*B_mag/(4.0*PI*rho) )+( abs(Bdir)*ce/sqrt(PI*rho) ))+
			  sqrt((ce*ce)+( B_mag*B_mag/(4.0*PI*rho) )-( abs(Bdir)*ce/sqrt(PI*rho) ))) + abs(udir);
		if (af > a_speed(0)){a_speed(0) = af;}
		}
	}
	PROTO_KERNEL_END(waveSpeedBoundF, waveSpeedBound)
  
  
	
	PROTO_KERNEL_START
	void dot_pro_calcFF(State& a_dot_pro,
					 const Var<double,1>& a_d_perp_N_s,
					 const State& a_d_perp_F)
	{
		for (int i=0; i< NUMCOMPS; i++){
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
		for (int i=0; i< NUMCOMPS; i++){
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
		static Stencil<double> m_ahead_shift[DIM];
		static Stencil<double> m_behind_shift[DIM];
		static Stencil<double> m_copy_f[DIM];
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
			m_ahead_shift[dir] = 1.0*Shift(Point::Basis(dir)*(1));
			m_behind_shift[dir] = 1.0*Shift(Point::Basis(dir)*(-1));
			m_copy_f[dir] = 1.0*Shift(Point::Zeros());
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

		
		Scalar Jacobian_ave(dbx0);
		MHD_Mapping::Jacobian_Ave_calc(Jacobian_ave,a_dx,dbx0);
		
		Vector a_U_demapped(dbx0);
		MHD_Mapping::JU_to_U_calc(a_U_demapped, a_U, Jacobian_ave, dbx0);
		
		
		//Vector W_bar(dbx0);
		//consToPrimcalc(W_bar,a_U_demapped,gamma);
		
		Vector W_bar = forall<double,NUMCOMPS>(consToPrim,a_U_demapped, gamma); 

		Vector U = m_deconvolve(a_U_demapped);

		Vector W    = forall<double,NUMCOMPS>(consToPrim,U, gamma);
		if(a_computeMaxWaveSpeed)
		{
		  Scalar umax = forall<double>(waveSpeedBound,a_rangeBox,W, gamma);
	 
		  umax.absMax(a_Rxn);  // returns 0 when used with CUDA
		}


		Vector W_ave = m_laplacian(W_bar,1.0/24.0);
		W_ave += W;
		
		bool limiter_apply = true;

		for (int d = 0; d < DIM; d++)
		{

			//Vector W_ave_low = m_interp_L[d](W_ave);
			//Vector W_ave_high = m_interp_H[d](W_ave);
		
			Vector W_ave_low = m_interp_edge[d](W_ave);
			Vector W_ave_high = m_interp_edge[d](W_ave);
			Vector W_ave_low_lim_flat(dbx0),W_ave_high_lim_flat(dbx0);
			if (limiter_apply){
				MHD_Limiters::MHD_Limiters(W_ave_low_lim_flat,W_ave_high_lim_flat,W_ave_low,W_ave_high,W_ave,d,a_dx);
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
			Scalar Lambda_f(dbx0), N_s_d_ave_f(dbx0);
			Vector F_f_mapped(dbx0);
			F_f_mapped.setVal(0.0);
			for (int s = 0; s < DIM; s++){
				Lambda_f = forall<double>(lambdacalc, W_low, W_high, s, gamma);
				F_low = forall<double,NUMCOMPS>(getFlux, W_low, s, gamma);
				F_high = forall<double,NUMCOMPS>(getFlux, W_high, s, gamma);
				F_f = forall<double,NUMCOMPS>(rusanovState, W_low, W_high, F_low, F_high, Lambda_f, s, gamma);
				MHD_Mapping::N_ave_f_calc_func(N_s_d_ave_f,s,d,a_dx);
#if DIM>1				
				F_ave_f = m_convolve_f[d](F_f);
#else
				F_ave_f = F_f;
#endif
				Vector dot_pro_sum(dbx0);
				dot_pro_sum.setVal(0.0);
				for (int s_temp = 0; s_temp < DIM; s_temp++){
					if (s_temp != d){
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
		
		}
		//PR_TIME("EulerOp::operator::RHS*=-1.0/dx");
		a_Rhs *= -1./a_dx;
    }

}
