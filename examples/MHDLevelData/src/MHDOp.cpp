#include "Proto.H"
#include "MHDOp.H"
#include "CommonTemplates.H"
#include "Proto_Timer.H"
#include "Proto_WriteBoxData.H"
// For Chrono Timer (Talwinder)
#include <chrono>
#include <iostream>

//////////////////////////////


typedef BoxData<double,1,1,1> Scalar;
typedef BoxData<double,NUMCOMPS,1,1> Vector;

namespace MHDOp {
	
	
	PROTO_KERNEL_START
	void 
	PowellF(State&         a_P, 
			const State&   a_W)
	{
#if DIM==1	
		a_P(0) = 0.;
		a_P(1) = a_W(3);
		a_P(2) = a_W(1);
		a_P(3) = a_W(1)*a_W(3);
#endif	

#if DIM==2
		a_P(0) = 0.;
		a_P(1) = a_W(4);
		a_P(2) = a_W(5);
		a_P(3) = a_W(1);
		a_P(4) = a_W(2);
		a_P(5) = a_W(1)*a_W(4) + a_W(2)*a_W(5);
#endif

#if DIM==3	
		a_P(0) = 0.;
		a_P(1) = a_W(5);
		a_P(2) = a_W(6);
		a_P(3) = a_W(7);
		a_P(4) = a_W(1);
		a_P(5) = a_W(2);
		a_P(6) = a_W(3);
		a_P(7) = a_W(1)*a_W(5) + a_W(2)*a_W(6) + a_W(3)*a_W(7);
#endif	
	}
	  PROTO_KERNEL_END(PowellF, Powell)
	
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
	void BavgcalcF(State& a_Bavg,
					 const State& a_W_ave,
					 int a_d)
	{
		for (int i=0; i< NUMCOMPS; i++){
			a_Bavg(i) = a_W_ave(2+DIM+a_d);
		}
	}
	PROTO_KERNEL_END(BavgcalcF, Bavgcalc) 


	PROTO_KERNEL_START
	void del_W_f_m_calcF(State& a_del_W_f_m,
					 const State& a_W_ave,
					 const State& a_W_ave_high)
	{
		for (int i=0; i< NUMCOMPS; i++){
			a_del_W_f_m(i) = a_W_ave(i) - a_W_ave_high(i);
		}
	}
	PROTO_KERNEL_END(del_W_f_m_calcF, del_W_f_m_calc) 

	PROTO_KERNEL_START
	void del_W_f_p_calcF(State& a_del_W_f_p,
					 const State& a_W_ave,
					 const State& a_W_ave_low_ahead)
	{
		for (int i=0; i< NUMCOMPS; i++){
			a_del_W_f_p(i) = a_W_ave_low_ahead(i) - a_W_ave(i);
		}
	}
	PROTO_KERNEL_END(del_W_f_p_calcF, del_W_f_p_calc)


	PROTO_KERNEL_START
	void del2_W_f_calcF(State& del2_W_f,
					 const State& a_W_ave,
					 const State& a_W_ave_high,
					 const State& a_W_ave_low_ahead)
	{
		for (int i=0; i< NUMCOMPS; i++){
			del2_W_f(i) = 6.0*(a_W_ave_high(i) - 2.0*a_W_ave(i) + a_W_ave_low_ahead(i));
		}
	}
	PROTO_KERNEL_END(del2_W_f_calcF, del2_W_f_calc)


	PROTO_KERNEL_START
	void del2_W_c_calcF(State& del2_W_c,
					 const State& a_W_ave,
					 const State& a_W_ave_behind,
					 const State& a_W_ave_ahead)
	{
		for (int i=0; i< NUMCOMPS; i++){
			del2_W_c(i) = a_W_ave_behind(i) - 2.0*a_W_ave(i) + a_W_ave_ahead(i);
		}
	}
	PROTO_KERNEL_END(del2_W_c_calcF, del2_W_c_calc)


	PROTO_KERNEL_START
	void del3_W_calcF(State& a_del3_W,
					 const State& a_del2_W_c,
					 const State& a_del2_W_c_behind)
	{
		for (int i=0; i< NUMCOMPS; i++){
			a_del3_W(i) = a_del2_W_c(i) - a_del2_W_c_behind(i);
		}
	}
	PROTO_KERNEL_END(del3_W_calcF, del3_W_calc)



	PROTO_KERNEL_START
	void del2_W_lim_calcF(State& a_del2_W_lim,
					 const State& a_del2_W_f,
					 const State& a_del2_W_c,
					 const State& a_del2_W_c_ahead,
					 const State& a_del2_W_c_behind)
	{
		for (int i=0; i< NUMCOMPS; i++){
			if ((a_del2_W_c_behind(i) >= 0.0) && (a_del2_W_c(i) >= 0.0) && (a_del2_W_c_ahead(i) >= 0.0) && (a_del2_W_f(i) >= 0.0)){
				a_del2_W_lim(i) = std::min({std::abs(a_del2_W_f(i)), 1.25*std::abs(a_del2_W_c_behind(i)),
				1.25*std::abs(a_del2_W_c(i)), 1.25*std::abs(a_del2_W_c_ahead(i))} );										
			} else if ((a_del2_W_c_behind(i) < 0.0) && (a_del2_W_c(i) < 0.0) && (a_del2_W_c_ahead(i) < 0.0) && (a_del2_W_f(i) < 0.0)){
				a_del2_W_lim(i) = -1.0*std::min({std::abs(a_del2_W_f(i)), 1.25*std::abs(a_del2_W_c_behind(i)),
				1.25*std::abs(a_del2_W_c(i)), 1.25*std::abs(a_del2_W_c_ahead(i))} );
			} else {
				a_del2_W_lim(i) = 0.0;
			}	
		}
	}
	PROTO_KERNEL_END(del2_W_lim_calcF, del2_W_lim_calc)


	PROTO_KERNEL_START
	void rho_i_calcF(State& a_rho_i,
					 const State& a_del2_W_f,
					 const State& a_del2_W_lim,
					 const State& a_W_ave,
					 const State& a_W_ave_ahead,
					 const State& a_W_ave_ahead2,
					 const State& a_W_ave_behind,
					 const State& a_W_ave_behind2)
	{
		for (int i=0; i< NUMCOMPS; i++){
			double rhs = 1.0e-12*std::max({std::abs(a_W_ave_behind2(i)), std::abs(a_W_ave_behind(i)), std::abs(a_W_ave(i)), std::abs(a_W_ave_ahead(i)), std::abs(a_W_ave_ahead2(i))});
			double lhs = std::abs(a_del2_W_f(i));
			if (lhs <= rhs){	
				a_rho_i(i) = 0.0;
			} else {
				a_rho_i(i) = a_del2_W_lim(i)/a_del2_W_f(i);
				//if (a_rho_i(i) < 0.99) cout << i << " " << a_del2_W_lim(i) << " " << a_del2_W_f(i) << " " << a_rho_i(i) << endl;
			}			
		}
	}
	PROTO_KERNEL_END(rho_i_calcF, rho_i_calc)


	PROTO_KERNEL_START
	void del3_W_min_calcF(State& a_del3_W_min,
					 const State& a_del3_W_behind,
					 const State& a_del3_W,
					 const State& a_del3_W_ahead,
					 const State& a_del3_W_ahead2)
	{
		for (int i=0; i< NUMCOMPS; i++){
			a_del3_W_min(i) = std::min({a_del3_W_behind(i), a_del3_W(i), a_del3_W_ahead(i), a_del3_W_ahead2(i)});			
		}
	}
	PROTO_KERNEL_END(del3_W_min_calcF, del3_W_min_calc)


	PROTO_KERNEL_START
	void del3_W_max_calcF(State& a_del3_W_max,
					 const State& a_del3_W_behind,
					 const State& a_del3_W,
					 const State& a_del3_W_ahead,
					 const State& a_del3_W_ahead2)
	{
		for (int i=0; i< NUMCOMPS; i++){
			a_del3_W_max(i) = std::max({a_del3_W_behind(i), a_del3_W(i), a_del3_W_ahead(i), a_del3_W_ahead2(i)});			
		}
	}
	PROTO_KERNEL_END(del3_W_max_calcF, del3_W_max_calc)

    const bool printlim = true;

	PROTO_KERNEL_START
	void limiter_low_calcF(
	                 const Point& a_pt,
					 State& a_W_ave_low_ahead_limited,
					 const State& a_W_ave_low_ahead,
					 const State& a_del_W_f_m,
					 const State& a_del_W_f_p,
					 const State& a_W_ave,
					 const State& a_W_ave_ahead2,
					 const State& a_W_ave_behind2,
					 const State& a_rho_i,
					 const State& a_del3_W_max,
					 const State& a_del3_W_min,
					 int   a_dir,
					 const double a_dx)
	{
		double rhs = 1.0 - 1.0e-12;
		double rhs_test, lhs_test, rhs_test2, lhs_test2;
		double x_loc = a_pt[0];//*a_dx + 0.5*a_dx;
		double y_loc = a_pt[1];//*a_dx + 0.5*a_dx;
		double x_loc2 = a_pt[0]*a_dx + 0.5*a_dx;
		double y_loc2 = a_pt[1]*a_dx + 0.5*a_dx;
		
		for (int i=0; i< NUMCOMPS; i++)
		{ 
	      if (i != (NUMCOMPS + a_dir - 2))
		  {	  
	        if ((a_del_W_f_m(i) * a_del_W_f_p(i) <= 0.0) || ((a_W_ave(i) - a_W_ave_behind2(i))*(a_W_ave_ahead2(i) - a_W_ave(i)) <= 0.0)){
				if (a_rho_i(i) < rhs){
					rhs_test = a_del3_W_max(i) - a_del3_W_min(i);
					lhs_test = 0.1*std::max({std::abs(a_del3_W_min(i)), std::abs(a_del3_W_max(i))});
					double prevval = a_W_ave_low_ahead(i);
					if (lhs_test <= rhs_test){
						
						rhs_test2 = 2.0*std::abs(a_del_W_f_m(i));
						lhs_test2 = std::abs(a_del_W_f_p(i));
						if (a_del_W_f_m(i) * a_del_W_f_p(i) < 0.0){
							a_W_ave_low_ahead_limited(i) = a_W_ave(i) + a_rho_i(i)*a_del_W_f_p(i);
							if (printlim){
								if (std::abs(a_W_ave_low_ahead_limited(i) - a_W_ave_low_ahead(i)) > 1.0e-14 && x_loc2 > 0.46 && x_loc2 < .54 && y_loc2 > 0.35 && y_loc2 < .41) cout << "x = " << x_loc << " y = " << y_loc << " Component = "<< i << " dir = " << a_dir <<" Condition : 1_a_L) Before limit: " << a_W_ave_low_ahead(i) << " Change: " << a_W_ave_low_ahead_limited(i) - a_W_ave_low_ahead(i) << " rho: " << a_rho_i(i) << endl;
								//if (std::abs(a_W_ave_low_ahead_limited(i) - a_W_ave_low_ahead(i)) > 1.0e-14 && x_loc2 > 0.46 && x_loc2 < .54 && y_loc2 > 0.35 && y_loc2 < .41) cout << lhs_test << " " << rhs_test << endl;
							}
						} else if (lhs_test2 >= rhs_test2){
							a_W_ave_low_ahead_limited(i) = a_W_ave(i) + 2.0*(1.0-a_rho_i(i))*a_del_W_f_m(i) + a_rho_i(i)*a_del_W_f_p(i);
							if (printlim){
								if (std::abs(a_W_ave_low_ahead_limited(i) - a_W_ave_low_ahead(i)) > 1.0e-14 && x_loc2 > 0.46 && x_loc2 < .54 && y_loc2 > 0.35 && y_loc2 < .41) cout << "x = " << x_loc << " y = " << y_loc << " Component = "<< i << " dir = " << a_dir <<" Condition : 1_c_L) Before limit: " << a_W_ave_low_ahead(i) << " Change: " << a_W_ave_low_ahead_limited(i) - a_W_ave_low_ahead(i) <<" rho: " << a_rho_i(i)<< endl;
								//if (std::abs(a_W_ave_low_ahead_limited(i) - a_W_ave_low_ahead(i)) > 1.0e-14 && x_loc2 > 0.46 && x_loc2 < .54 && y_loc2 > 0.35 && y_loc2 < .41) cout << lhs_test << " " << rhs_test << endl;
							}
						}
					}
				}
			} else {
				rhs_test2 = 2.0*std::abs(a_del_W_f_m(i));
				lhs_test2 = std::abs(a_del_W_f_p(i));
				if (lhs_test2 >= rhs_test2){
					a_W_ave_low_ahead_limited(i) = a_W_ave(i) + 2.0*a_del_W_f_m(i);
					if (printlim){
						if (std::abs(a_W_ave_low_ahead_limited(i) - a_W_ave_low_ahead(i)) > 1.0e-14 && x_loc2 > 0.46 && x_loc2 < .54 && y_loc2 > 0.35 && y_loc2 < .41) cout << "x = " << x_loc << " y = " << y_loc << " Component = "<< i << " dir = " << a_dir <<" Condition : 2_b_L) Before limit: " << a_W_ave_low_ahead(i) << " Change: " << a_W_ave_low_ahead_limited(i) - a_W_ave_low_ahead(i) << " rho: " << a_rho_i(i) << endl;
					}
				}
			}		
		  } 
		}
	}
	PROTO_KERNEL_END(limiter_low_calcF, limiter_low_calc)


	PROTO_KERNEL_START
	void limiter_high_calcF(
	                 const Point& a_pt, 
	                 State& a_W_ave_high_limited,
					 const State& a_W_ave_high,
					 const State& a_del_W_f_m,
					 const State& a_del_W_f_p,
					 const State& a_W_ave,
					 const State& a_W_ave_ahead2,
					 const State& a_W_ave_behind2,
					 const State& a_rho_i,
					 const State& a_del3_W_max,
					 const State& a_del3_W_min,
					 int   a_dir,
					 const double a_dx)
	{
		double rhs = 1.0 - 1.0e-12;
		double rhs_test, lhs_test, rhs_test2, lhs_test2;
		double x_loc = a_pt[0];//*a_dx + 0.5*a_dx;
		double y_loc = a_pt[1];//*a_dx + 0.5*a_dx;
		double x_loc2 = a_pt[0]*a_dx + 0.5*a_dx;
		double y_loc2 = a_pt[1]*a_dx + 0.5*a_dx;
		
		
		for (int i=0; i< NUMCOMPS; i++)
		{
		  if (i != (NUMCOMPS + a_dir - 2))
		  {	
	        if ((a_del_W_f_m(i) * a_del_W_f_p(i) <= 0.0) || ((a_W_ave(i) - a_W_ave_behind2(i))*(a_W_ave_ahead2(i) - a_W_ave(i)) <= 0.0)){
				if (a_rho_i(i) < rhs){
					rhs_test = a_del3_W_max(i) - a_del3_W_min(i);
					lhs_test = 0.1*std::max({std::abs(a_del3_W_min(i)), std::abs(a_del3_W_max(i))});
					if (lhs_test <= rhs_test){
						
						rhs_test2 = 2.0*std::abs(a_del_W_f_p(i));
						lhs_test2 = std::abs(a_del_W_f_m(i));
						if (a_del_W_f_m(i) * a_del_W_f_p(i) < 0.0){
							a_W_ave_high_limited(i) = a_W_ave(i) - a_rho_i(i)*a_del_W_f_m(i);
							if (printlim){
								if (std::abs(a_W_ave_high_limited(i) - a_W_ave_high(i)) > 1.0e-14 && x_loc2 > 0.46 && x_loc2 < .54 && y_loc2 > 0.35 && y_loc2 < .41) cout << "x = " << x_loc << " y = " << y_loc << " Component = "<< i << " dir = " << a_dir <<" Condition : 1_a_R) Before limit: " << a_W_ave_high(i) << " Change: " << a_W_ave_high_limited(i) - a_W_ave_high(i) << " rho: " << a_rho_i(i) << endl;
								//if (std::abs(a_W_ave_high_limited(i) - a_W_ave_high(i)) > 1.0e-14 && x_loc2 > 0.46 && x_loc2 < .54 && y_loc2 > 0.35 && y_loc2 < .41) cout << lhs_test << " " << rhs_test << endl;
							}
						} else if (lhs_test2 >= rhs_test2){
							a_W_ave_high_limited(i) = a_W_ave(i) - 2.0*(1.0-a_rho_i(i))*a_del_W_f_p(i) - a_rho_i(i)*a_del_W_f_m(i);
							if (printlim){
								if (std::abs(a_W_ave_high_limited(i) - a_W_ave_high(i)) > 1.0e-14 && x_loc2 > 0.46 && x_loc2 < .54 && y_loc2 > 0.35 && y_loc2 < .41) cout << "x = " << x_loc << " y = " << y_loc << " Component = "<< i << " dir = " << a_dir <<" Condition : 1_b_R) Before limit: " << a_W_ave_high(i) << " Change: " << a_W_ave_high_limited(i) - a_W_ave_high(i) << " rho: " << a_rho_i(i) << " " << lhs_test2 << " " << rhs_test2/2.0 << endl;
								//if (std::abs(a_W_ave_high_limited(i) - a_W_ave_high(i)) > 1.0e-14 && x_loc2 > 0.46 && x_loc2 < .54 && y_loc2 > 0.35 && y_loc2 < .41) cout << lhs_test << " " << rhs_test << endl;
							}
						}
					}
				}
			} else {
				rhs_test2 = 2.0*std::abs(a_del_W_f_p(i));
				lhs_test2 = std::abs(a_del_W_f_m(i));
				if (lhs_test2 >= rhs_test2){
					a_W_ave_high_limited(i) = a_W_ave(i) - 2.0*a_del_W_f_p(i);
					if (printlim){
						if (std::abs(a_W_ave_high_limited(i) - a_W_ave_high(i)) > 1.0e-14 && x_loc2 > 0.46 && x_loc2 < .54 && y_loc2 > 0.35 && y_loc2 < .41) cout << "x = " << x_loc << " y = " << y_loc << " Component = "<< i << " dir = " << a_dir <<" Condition : 2_a_R) Before limit: " << a_W_ave_high(i) << " Change: " << a_W_ave_high_limited(i) - a_W_ave_high(i) << " rho: " << a_rho_i(i) << endl;
					}
				}
			}		
		  }
		}
	}
	PROTO_KERNEL_END(limiter_high_calcF, limiter_high_calc)



	PROTO_KERNEL_START
	void eta_tilde_d_calcF(Var<double,1>& a_eta_tilde_d,
					 const State& a_W_ave,
					 const State& a_W_ave_ahead,
					 const State& a_W_ave_ahead2,
					 const State& a_W_ave_behind,
					 const State& a_W_ave_behind2,
					 int a_d)
	{   
		double z1=0.85, z0=0.75, delta = 0.33, si;
		double p, p_ahead, p_ahead2, p_behind, p_behind2;
#if DIM==1	
		p = a_W_ave(2) + a_W_ave(3)*a_W_ave(3)/8.0/PI;
		p_ahead = a_W_ave_ahead(2) + a_W_ave_ahead(3)*a_W_ave_ahead(3)/8.0/PI;
		p_ahead2 = a_W_ave_ahead2(2) + a_W_ave_ahead2(3)*a_W_ave_ahead2(3)/8.0/PI;
		p_behind = a_W_ave_behind(2) + a_W_ave_behind(3)*a_W_ave_behind(3)/8.0/PI;
		p_behind2 = a_W_ave_behind2(2) + a_W_ave_behind2(3)*a_W_ave_behind2(3)/8.0/PI;
#endif
#if DIM==2	
		p = a_W_ave(3) + (a_W_ave(4)*a_W_ave(4)+a_W_ave(5)*a_W_ave(5))/8.0/PI;
		p_ahead = a_W_ave_ahead(3) + (a_W_ave_ahead(4)*a_W_ave_ahead(4)+a_W_ave_ahead(5)*a_W_ave_ahead(5))/8.0/PI;
		p_ahead2 = a_W_ave_ahead2(3) + (a_W_ave_ahead2(4)*a_W_ave_ahead2(4)+a_W_ave_ahead2(5)*a_W_ave_ahead2(5))/8.0/PI;
		p_behind = a_W_ave_behind(3) + (a_W_ave_behind(4)*a_W_ave_behind(4)+a_W_ave_behind(5)*a_W_ave_behind(5))/8.0/PI;
		p_behind2 = a_W_ave_behind2(3) + (a_W_ave_behind2(4)*a_W_ave_behind2(4)+a_W_ave_behind2(5)*a_W_ave_behind2(5))/8.0/PI;
#endif
#if DIM==3
		p = a_W_ave(4) + (a_W_ave(5)*a_W_ave(5)+a_W_ave(6)*a_W_ave(6)+a_W_ave(7)*a_W_ave(7))/8.0/PI;
		p_ahead = a_W_ave_ahead(4) + (a_W_ave_ahead(5)*a_W_ave_ahead(5)+a_W_ave_ahead(6)*a_W_ave_ahead(6)+a_W_ave_ahead(7)*a_W_ave_ahead(7))/8.0/PI;
		p_ahead2 = a_W_ave_ahead2(4) + (a_W_ave_ahead2(5)*a_W_ave_ahead2(5)+a_W_ave_ahead2(6)*a_W_ave_ahead2(6)+a_W_ave_ahead2(7)*a_W_ave_ahead2(7))/8.0/PI;
		p_behind = a_W_ave_behind(4) + (a_W_ave_behind(5)*a_W_ave_behind(5)+a_W_ave_behind(6)*a_W_ave_behind(6)+a_W_ave_behind(7)*a_W_ave_behind(7))/8.0/PI;
		p_behind2 = a_W_ave_behind2(4) + (a_W_ave_behind2(5)*a_W_ave_behind2(5)+a_W_ave_behind2(6)*a_W_ave_behind2(6)+a_W_ave_behind2(7)*a_W_ave_behind2(7))/8.0/PI;
#endif
		double arg = std::abs(p_ahead-p_behind)/std::abs(p_ahead2-p_behind2);
		if (arg > z1){si = 0.0;} 
		if (arg < z0){si = 1.0;} 
		if ((arg<=z1) && (arg>=z0)){si = 1.0 - (arg-z0)/(z1-z0);}
		
		double lhs1 = a_W_ave_behind(1+a_d)-a_W_ave_ahead(1+a_d);
		double lhs2 = std::abs(p_ahead-p_behind)/std::min({p_ahead,p_behind});
		if ((lhs1 > 0.0) && (lhs2 > delta)){
			a_eta_tilde_d(0) = si;
		} else {
			a_eta_tilde_d(0) = 0.0;
		}
	}
	PROTO_KERNEL_END(eta_tilde_d_calcF, eta_tilde_d_calc) 

	PROTO_KERNEL_START
	void eta_d_calcF(Var<double,1>& a_eta_d,
					 const Var<double,1>& a_eta_tilde_d,
					 const Var<double,1>& a_eta_tilde_d_ahead,
					 const Var<double,1>& a_eta_tilde_d_behind)
	{   

		a_eta_d(0) = std::min({a_eta_tilde_d(0),a_eta_tilde_d_ahead(0),a_eta_tilde_d_behind(0)});

	}
	PROTO_KERNEL_END(eta_d_calcF, eta_d_calc) 

	PROTO_KERNEL_START
	void eta_calcF(Var<double,1>& a_eta,
					 const Var<double,1>& a_eta_d,
					 const Var<double,1>& a_eta_d_old)
	{   

		a_eta(0) = std::min({a_eta_d(0),a_eta_d_old(0)});

	}
	PROTO_KERNEL_END(eta_calcF, eta_calc) 

	PROTO_KERNEL_START
	void Flat_calcF(State& a_flattened,
					 const State& a_not_flattened,
					 const State& a_W_ave,
					 const Var<double,1>& a_eta,
					 int   a_dir)
	{   
		
		for (int i=0; i< NUMCOMPS; i++){
			if (i != (NUMCOMPS + a_dir - 2))
	        //if (i != 100)
		    //if (i < 4)
		    {
				a_flattened(i) = (1.0-a_eta(0))*a_not_flattened(i) + a_eta(0)*a_W_ave(i);	
		    } else {
				a_flattened(i) = a_not_flattened(i);
			}				
		}

	}
	PROTO_KERNEL_END(Flat_calcF, Flat_calc) 

	PROTO_KERNEL_START
	void v_d_calcF(Var<double,1>& a_v_d,
					 const State& a_W_bar,
					 const int a_d)
	{
			a_v_d(0) = a_W_bar(1+a_d);			
	}
	PROTO_KERNEL_END(v_d_calcF, v_d_calc) 

	PROTO_KERNEL_START
	void viscosity1_calcF(Var<double,1>& a_viscosity,
					 const Var<double,1>& a_v,
					 const Var<double,1>& a_v_behind)
	{
			a_viscosity(0) = a_v(0)-a_v_behind(0);			
	}
	PROTO_KERNEL_END(viscosity1_calcF, viscosity1_calc) 

	PROTO_KERNEL_START
	void v_d2_div_calcF(Var<double,1>& a_v_d2_div,
					 const Var<double,1>& v_d2_ahead,
					 const Var<double,1>& v_d2_behind,
					 const Var<double,1>& v_d2_behind_dp,
					 const Var<double,1>& v_d2_behind_dm)
	{
			a_v_d2_div(0) = (v_d2_ahead(0)-v_d2_behind(0)+v_d2_behind_dp(0)-v_d2_behind_dm(0))/4.0;			
	}
	PROTO_KERNEL_END(v_d2_div_calcF, v_d2_div_calc) 


	PROTO_KERNEL_START
	void Fast_MS_speed_calcF(Var<double,1>& a_Fast_MS_speed,
					 const State& a_W_bar,
					 int a_d,
					 double a_gamma)
	{
		double gamma = a_gamma;
		double rho=0., p=0., Bx=0., By=0., Bz=0., ce, B_mag, Bdir;
#if DIM == 1
		rho = a_W_bar(0);
		p   = a_W_bar(2);
		Bx  = a_W_bar(3);
#endif	
#if DIM == 2
		rho = a_W_bar(0);
		p   = a_W_bar(3);
		Bx  = a_W_bar(4);
		By  = a_W_bar(5);
#endif	
#if DIM == 3
		rho = a_W_bar(0);
		p   = a_W_bar(4);
		Bx  = a_W_bar(5);
		By  = a_W_bar(6);
		Bz  = a_W_bar(7);
#endif	
		if (a_d == 0){
			Bdir = Bx;
		};
		if (a_d == 1){
			Bdir = By;
		};
		if (a_d == 2){
			Bdir = Bz;
		};
		
		ce = sqrt(gamma*p/rho);
		B_mag = sqrt(Bx*Bx+By*By+Bz*Bz);
		a_Fast_MS_speed(0) = 0.5*(sqrt((ce*ce)+( B_mag*B_mag/(4.0*PI*rho) )+( abs(Bdir)*ce/sqrt(PI*rho) ))+
			  sqrt((ce*ce)+( B_mag*B_mag/(4.0*PI*rho) )-( abs(Bdir)*ce/sqrt(PI*rho) )));
	}
	PROTO_KERNEL_END(Fast_MS_speed_calcF, Fast_MS_speed_calc)

	PROTO_KERNEL_START
	void Fast_MS_speed_min_calcF(Var<double,1>& a_Fast_MS_speed_min,
					 const Var<double,1>& a_Fast_MS_speed,
					 const Var<double,1>& a_Fast_MS_speed_behind)
	{   
		a_Fast_MS_speed_min(0) = std::min({a_Fast_MS_speed(0),a_Fast_MS_speed_behind(0)});
	}
	PROTO_KERNEL_END(Fast_MS_speed_min_calcF, Fast_MS_speed_min_calc) 


	PROTO_KERNEL_START
	void Visc_coef_calcF(Var<double,1>& a_Visc_coef,
					 const Var<double,1>& a_h_lambda,
					 const Var<double,1>& a_Fast_MS_speed_min)
	{   
		if (a_h_lambda(0) <= 0){
			double temp = a_h_lambda(0)*a_h_lambda(0)/0.3/a_Fast_MS_speed_min(0)/a_Fast_MS_speed_min(0);
			a_Visc_coef(0) = a_h_lambda(0)*std::min({temp,1.0});
		} else {
			a_Visc_coef(0) = 0.0;
		}
	}
	PROTO_KERNEL_END(Visc_coef_calcF, Visc_coef_calc) 

	PROTO_KERNEL_START
	void mu_f_calcF(State& a_mu_f,
					 const Var<double,1>& a_Visc_coef,
					 const State& a_U,
					 const State& a_U_behind)
	{   
		for (int i=0; i< NUMCOMPS; i++){
			a_mu_f(i) = 0.3*a_Visc_coef(0)*(a_U(i)-a_U_behind(i));
		}
	}
	PROTO_KERNEL_END(mu_f_calcF, mu_f_calc) 


 

PROTO_KERNEL_START 
	unsigned int InitializeStateF(State& a_U,
                             const V& a_x,
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
	
	
	// //////Modifying parameters for 2D current sheet problem/////
	// rho = 1.0;
    // p = 0.1;
	// u = 0.1 * sin(2*PI*a_x(1));
	// if (a_x(0) >= 0.0 && a_x(0) < 0.5){By = 1.0;}
	// if (a_x(0) >= 0.5 && a_x(0) < 1.5){By = -1.0;}
	// if (a_x(0) >= 1.5 && a_x(0) <= 2.0){By = 1.0;}
	
	// rho = 1.0;
    // p = 0.1;
	// v = 0.1 * sin(2.0*PI*a_x(0));
	// if (a_x(1) >= 0.0 && a_x(1) < 0.5){Bx = 1.0;}
	// if (a_x(1) >= 0.5 && a_x(1) < 1.5){Bx = -1.0;}
	// if (a_x(1) >= 1.5 && a_x(1) <= 2.0){Bx = 1.0;}
	// ////////////////////////////////////////////////////////////
	
	
	//////Modifying parameters for 2D Orszag Tang problem///////
	// Case 1:
	rho = gamma*((2.0 * gamma)/8.0/PI)*1.0;
    p = (2.0 * gamma)/8.0/M_PI;
	u = -sin(2.0*PI*a_x(1));
	v =  sin(2.0*PI*a_x(0));
	Bx = -sin(2.0*PI*a_x(1));
	By =  sin(4.0*PI*a_x(0));
	// Case 2:
	// rho = 1.0;
    // p = 1.0/gamma;
	// u = -sin(2.0*PI*a_x(1));
	// v =  sin(2.0*PI*a_x(0));
	// Bx = -sin(2.0*PI*a_x(1))/gamma;
	// By =  sin(4.0*PI*a_x(0))/gamma;
	////////////////////////////////////////////////////////////
	
	
	//////Modifying parameters for Alfven wave problem///////
	// rho = 1.0;
    // p = 1.0;
	// u =  sin(2.0*PI*a_x(0));
	// Bx = sin(2.0*PI*a_x(0));
	////////////////////////////////////////////////////////////

    //////Modifying parameters for Acoustic pulse problem///////
	// double rho_0 = 1.4;
	// double delta_rho_0 = 0.14;
	// double rad = sqrt((a_x(0)-0.5)*(a_x(0)-0.5) + (a_x(1)-0.5)*(a_x(1)-0.5));
	// if (rad < 0.5){
		// rho = rho_0 + delta_rho_0*exp(-16.0*rad*rad)*pow(cos(PI*rad),6.0);
	// } else {
		// rho = rho_0;
	// }
	// p = pow(rho/rho_0,gamma);
	////////////////////////////////////////////////////////////	
	
	
	//////Modifying parameters for Euler problem///////
	// double rho0 = gamma;
    // double p0 = 1.;
    // rho = rho0;
    // rho += .01*rho0*sin(2*2*PI*a_x(0));
    // p = p0*pow(rho/rho0,gamma);
    // double c0 = sqrt(gamma*p0/rho0);
    // double c = sqrt(gamma*p/rho);
    // u = 2*(c-c0)/(gamma-1.);
	////////////////////////////////////////////////////////////
	
	
	
	
	
	double e = p/(gamma-1.0) + rho*(u*u+v*v+w*w)/2.0 + (Bx*Bx+By*By+Bz*Bz)/8.0/PI;
	
#if DIM == 1	
	a_U(0) = rho;      //rho
	a_U(1) = rho*u;    //Momentum-x
	a_U(2) = e;        //Energy
	a_U(3) = Bx;	   //Bx
#endif	
#if DIM == 2
	a_U(0) = rho;      //rho
	a_U(1) = rho*u;    //Momentum-x
	a_U(2) = rho*v;    //Momentum-y 
	a_U(3) = e;        //Energy
	a_U(4) = Bx;	   //Bx
	a_U(5) = By;       //By
#endif
#if DIM == 3
	a_U(0) = rho;      //rho
	a_U(1) = rho*u;    //Momentum-x
	a_U(2) = rho*v;    //Momentum-y 
	a_U(3) = rho*w;    //Momentum-z 
	a_U(4) = e;        //Energy
	a_U(5) = Bx;	   //Bx
	a_U(6) = By;       //By
	a_U(7) = Bw;       //Bz
#endif
	
    return 0;
	}
	PROTO_KERNEL_END(InitializeStateF, InitializeState)




PROTO_KERNEL_START
void iotaFuncF(Point           & a_p,
               V               & a_X,
               double            a_h)
{
  for (int ii = 0; ii < DIM; ii++)
  {
    a_X(ii) = a_p[ii]*a_h + 0.5*a_h;
  }
}
PROTO_KERNEL_END(iotaFuncF,iotaFunc)

void initializeState(BoxData<double,NUMCOMPS>& a_state,
                     const double a_dx,
                     const double a_gamma)
{
    Box dbx0=a_state.box();
    Box dbx = dbx0.grow(NGHOST);
    Box dbx1 = dbx.grow(1);
    BoxData<double,NUMCOMPS> UBig(dbx1);
    BoxData<double,DIM> x(dbx1);
    forallInPlace_p(iotaFunc, dbx1, x, a_dx);
    forallInPlace(InitializeState,dbx1,UBig,x,a_gamma);
    //cout << "1D in 1: " << test1D(UBig,1) << endl;
    Stencil<double> Lap2nd = Stencil<double>::Laplacian();
    a_state |= Lap2nd(UBig,dbx,1.0/24.0); 
    a_state += UBig;
    //cout << "1D in 1: " << test1D(a_state,1) << endl;
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
    static Stencil<double> m_laplacian;
		static Stencil<double> m_deconvolve;
		static Stencil<double> m_copy;
		static Stencil<double> m_laplacian_f[DIM];
		static Stencil<double> m_deconvolve_f[DIM];
		static Stencil<double> m_convolve_f[DIM];
		static Stencil<double> m_interp_H[DIM];
		static Stencil<double> m_interp_L[DIM];
		static Stencil<double> m_interp_edge[DIM];
		//static Stencil<double> m_invert[DIM];
		static Stencil<double> m_divergence[DIM];
		static Stencil<double> m_ahead_shift[DIM];
		static Stencil<double> m_behind_shift[DIM];
		static Stencil<double> m_copy_f[DIM];
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
			//m_invert[dir] = Stencil<double>::invert(dir);
			m_divergence[dir] = Stencil<double>::FluxDivergence(dir);
			m_ahead_shift[dir] = 1.0*Shift(Point::Basis(dir)*(1));
			m_behind_shift[dir] = 1.0*Shift(Point::Basis(dir)*(-1));
			m_copy_f[dir] = 1.0*Shift(Point::Zeros());
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
    Vector W_bar = forall<double,NUMCOMPS>(consToPrim,a_U, gamma); // 

    Vector U = m_deconvolve(a_U);

    Vector W    = forall<double,NUMCOMPS>(consToPrim,U, gamma);
    if(a_computeMaxWaveSpeed)
    {
      Scalar umax = forall<double>(waveSpeedBound,a_rangeBox,W, gamma);
 
      umax.absMax(a_Rxn);  // returns 0 when used with CUDA
    }


    Vector W_ave = m_laplacian(W_bar,1.0/24.0);
    W_ave += W;
	
	bool limiter_apply = true;
	Scalar eta ;
		Scalar eta_old ;
	if (limiter_apply){
	
		
		
		
		for (int d = 0; d < DIM; d++)
		{

			Vector W_ave_ahead = alias(W_ave,Point::Basis(d)*(-1));
			Vector W_ave_ahead2 = alias(W_ave,Point::Basis(d)*(-2));
			Vector W_ave_ahead3 = alias(W_ave,Point::Basis(d)*(-3));
			Vector W_ave_behind = alias(W_ave, Point::Basis(d)*(1));
			Vector W_ave_behind2 = alias(W_ave,Point::Basis(d)*(2));	
			Vector W_ave_behind3 = alias(W_ave,Point::Basis(d)*(3));
			
			Scalar eta_tilde_d = forall<double>( eta_tilde_d_calc, W_ave, W_ave_ahead, W_ave_ahead2,
										  W_ave_behind, W_ave_behind2, d);
									  
			Scalar eta_tilde_d_ahead = alias(eta_tilde_d,Point::Basis(d)*(-1));
			Scalar eta_tilde_d_behind = alias(eta_tilde_d,Point::Basis(d)*(1));	
			
			Scalar eta_d = forall<double>( eta_d_calc, eta_tilde_d, eta_tilde_d_ahead, eta_tilde_d_behind);
			if (d>0){
				eta = forall<double>(eta_calc, eta_d, eta_old);
			} else {
				eta = m_copy(eta_d);
			}	
			eta_old = m_copy(eta);
		}
	}
	
    for (int d = 0; d < DIM; d++)
      {

        //Vector W_ave_low = m_interp_L[d](W_ave);
        //Vector W_ave_high = m_interp_H[d](W_ave);
		
		Vector W_ave_low = m_interp_edge[d](W_ave);
        Vector W_ave_high = m_interp_edge[d](W_ave);
		 
		 
		 
		 
		 if (limiter_apply){
		 
		 //Limiter Starts here
			Vector W_ave_low_ahead = alias(W_ave_low,Point::Basis(d)*(-1));
			Vector del_W_f_m = forall<double,NUMCOMPS>( del_W_f_m_calc, W_ave, W_ave_high);
			Vector del_W_f_p = forall<double,NUMCOMPS>(del_W_f_p_calc, W_ave, W_ave_low_ahead);
			Vector del2_W_f  = forall<double,NUMCOMPS>(del2_W_f_calc, W_ave, W_ave_high, W_ave_low_ahead);
			Vector W_ave_ahead = alias(W_ave,Point::Basis(d)*(-1));
			Vector W_ave_ahead2 = alias(W_ave,Point::Basis(d)*(-2));
			Vector W_ave_behind = alias(W_ave,Point::Basis(d)*(1));
			Vector W_ave_behind2 = alias(W_ave,Point::Basis(d)*(2));
			Vector del2_W_c  = forall<double,NUMCOMPS>( del2_W_c_calc, W_ave, W_ave_behind, W_ave_ahead);
			
			
			
			Vector del2_W_c_ahead = alias(del2_W_c,Point::Basis(d)*(-1));
			Vector del2_W_c_behind = alias(del2_W_c,Point::Basis(d)*(1));
			Vector del3_W = forall<double,NUMCOMPS>( del3_W_calc, del2_W_c, del2_W_c_behind);
			
			// if (d==0){
				// std::string filename="del3Vx_dx3_";
				// BoxData<double> del3_Vx = slice(del3_W,1);
				// WriteBoxData(filename.c_str(),del3_Vx,a_dx);
			// }
			if (d==0){
				std::string filename="W_ave_high";
				BoxData<double> W_ave_high_out = slice(W_ave_high,1);
				WriteBoxData(filename.c_str(),W_ave_high_out,a_dx);
			}
			if (d==0){
				std::string filename="W_ave";
				BoxData<double> W_ave_out = slice(W_ave,1);
				WriteBoxData(filename.c_str(),W_ave_out,a_dx);
			}
			
			Vector del3_W_behind = alias(del3_W,Point::Basis(d)*(1));
			Vector del3_W_ahead = alias(del3_W,Point::Basis(d)*(-1));
			Vector del3_W_ahead2 = alias(del3_W,Point::Basis(d)*(-2));
			Vector del2_W_lim = forall<double,NUMCOMPS>(del2_W_lim_calc, del2_W_f, del2_W_c, del2_W_c_ahead, del2_W_c_behind);
			Vector rho_i = forall<double,NUMCOMPS>( rho_i_calc, del2_W_f, del2_W_lim, W_ave, 
													W_ave_ahead, W_ave_ahead2, W_ave_behind, W_ave_behind2);
			Vector del3_W_min = forall<double,NUMCOMPS>(del3_W_min_calc, del3_W_behind, 
														 del3_W, del3_W_ahead, del3_W_ahead2);	
			Vector del3_W_max = forall<double,NUMCOMPS>( del3_W_max_calc, del3_W_behind, 
														 del3_W, del3_W_ahead, del3_W_ahead2);				 
			// Vector W_ave_low_ahead_limited = forall<double,NUMCOMPS>( limiter_low_calc, W_ave_low_ahead, del_W_f_m, del_W_f_p, W_ave, W_ave_ahead2, W_ave_behind2, rho_i, del3_W_max, del3_W_min, d);															  		
			// Vector W_ave_high_limited = forall<double,NUMCOMPS>( limiter_high_calc, W_ave_high, del_W_f_m, del_W_f_p, W_ave, W_ave_ahead2, W_ave_behind2, rho_i, del3_W_max, del3_W_min, d);	
			
			Vector W_ave_low_ahead_limited = m_copy_f[d](W_ave_low_ahead);
			Vector W_ave_high_limited = m_copy_f[d](W_ave_high);
            forallInPlace_p( limiter_low_calc,W_ave_low_ahead_limited, W_ave_low_ahead, del_W_f_m, del_W_f_p, W_ave, W_ave_ahead2, W_ave_behind2, rho_i, del3_W_max, del3_W_min, d, a_dx);
			forallInPlace_p( limiter_high_calc, W_ave_high_limited, W_ave_high, del_W_f_m, del_W_f_p, W_ave, W_ave_ahead2, W_ave_behind2, rho_i, del3_W_max, del3_W_min, d, a_dx);			
			// Slope flattening starts here
			Vector W_ave_low_ahead_lim_flat = forall<double,NUMCOMPS>(Flat_calc, W_ave_low_ahead_limited, W_ave, eta, d);
			Vector W_ave_high_lim_flat = forall<double,NUMCOMPS>(Flat_calc, W_ave_high_limited, W_ave, eta, d);
			//Slope flattening ends here

			W_ave_low = alias(W_ave_low_ahead_lim_flat,Point::Basis(d)*(1));
			W_ave_high = m_copy_f[d](W_ave_high_lim_flat);
			
			// W_ave_low = alias(W_ave_low_ahead_limited,Point::Basis(d)*(1));
			// W_ave_high = m_copy_f[d](W_ave_high_limited);
			//Limiter ends here
		 
		 }
		 
		 
		 
		 
		 
		 
				
		
#if DIM>1
			Vector W_low = m_deconvolve_f[d](W_ave_low);		
			Vector W_high = m_deconvolve_f[d](W_ave_high);
#else
			Vector W_low = W_ave_low;
			Vector W_high = W_ave_high;
#endif	

			Scalar Lambda_f = forall<double>(lambdacalc, W_low, W_high, d, gamma);
			Vector F_low = forall<double,NUMCOMPS>(getFlux, W_low, d, gamma);
			Vector F_high = forall<double,NUMCOMPS>(getFlux, W_high, d, gamma);
			Vector F_f = forall<double,NUMCOMPS>(rusanovState, W_low, W_high, F_low, F_high, Lambda_f, d, gamma);
#if DIM>1				
			Vector F_ave_f = m_convolve_f[d](F_f);
#else
			Vector F_ave_f = F_f;
#endif		
		
        a_Rhs += m_divergence[d](F_ave_f);
	
      }
    //PR_TIME("EulerOp::operator::RHS*=-1.0/dx");
    a_Rhs *= -1./a_dx;
  }
  
  // Used to find the Powell term
	void step2(BoxData<double,NUMCOMPS>& a_Rhs,
            const BoxData<double,NUMCOMPS>& a_U,
            const Box& a_rangeBox,
            const double a_dx,
            const double a_gamma,
            Reduction<double>& a_Rxn,
            bool a_computeMaxWaveSpeed,
            bool a_callBCs)
	{
		static Stencil<double> m_laplacian;
		static Stencil<double> m_deconvolve;
		static Stencil<double> m_divergence[DIM];
		static bool initialized = false;
		if(!initialized)
		{
		  m_laplacian = Stencil<double>::Laplacian();
		  m_deconvolve = (-1.0/24.0)*m_laplacian + (1.0)*Shift(Point::Zeros());
		  for (int dir = 0; dir < DIM; dir++)
		  {
			m_divergence[dir] = Stencil<double>::FluxDivergence(dir);
		  }
		  initialized =  true;
		}


		using namespace std;
		//PR_TIME("MHDOp::operator");
		a_Rhs.setVal(0.0);
		Vector divB_term;
		divB_term *= 0.0;
		double gamma = a_gamma;
		double retval;

		Vector W_bar = forall<double,NUMCOMPS>(consToPrim,a_U, gamma);
	    Vector U = m_deconvolve(a_U);
	    Vector W = forall<double,NUMCOMPS>( consToPrim,U, gamma);
	    Vector W_ave = m_laplacian(W_bar);
	    W_ave *= (1.0/24.0);
	    W_ave += W;
	    Vector Powell_term = forall<double,NUMCOMPS>( Powell,W_ave);
	    for (int d = 0; d < DIM; d++)
	    {
			Vector B_ave = forall<double,NUMCOMPS>( Bavgcalc, W_ave, d);
			divB_term += m_divergence[d](B_ave);
			divB_term *= Powell_term;
			a_Rhs += divB_term;
	    }
		a_Rhs *= -1./a_dx;
	}
	
	
	
	// Used to implement artificial viscosity
	void step3(BoxData<double,NUMCOMPS>& a_Rhs,
            const BoxData<double,NUMCOMPS>& a_U,
            const Box& a_rangeBox,
            const double a_dx,
            const double a_gamma,
            Reduction<double>& a_Rxn,
            bool a_computeMaxWaveSpeed,
            bool a_callBCs)
	{

		static Stencil<double> m_divergence[DIM];
		static Stencil<double> m_behind_shift[DIM];
		static bool initialized = false;
		if(!initialized)
		{
		  for (int dir = 0; dir < DIM; dir++)
		  {
			m_divergence[dir] = Stencil<double>::FluxDivergence(dir);
		    m_behind_shift[dir] = 1.0*Shift(Point::Basis(dir)*(-1));
		  }
		  initialized =  true;
		}


		using namespace std;
		//PR_TIME("MHDOp::operator");
		a_Rhs.setVal(0.0);
		double gamma = a_gamma;
		double retval;

		Vector W_bar = forall<double,NUMCOMPS>( consToPrim,a_U, gamma);
	    for (int d = 0; d < DIM; d++)
	    { 
			Scalar v_d =  forall<double>(v_d_calc,W_bar,d);
			Scalar v_d_behind = alias(v_d,Point::Basis(d)*(1));
			Scalar h_lambda = forall<double>(viscosity1_calc,v_d,v_d_behind); 
			for (int d2 = 0; d2 < DIM; d2++){	
				if (d2!=d){
					Scalar v_d2 = forall<double>(v_d_calc,W_bar,d2); 
					Scalar v_d2_ahead = alias(v_d2,Point::Basis(d2)*(-1));
					Scalar v_d2_behind = alias(v_d2,Point::Basis(d2)*(1));
					Scalar v_d2_behind_dp = alias(v_d2_ahead,Point::Basis(d)*(1));
					Scalar v_d2_behind_dm = alias(v_d2_behind,Point::Basis(d)*(1));
					Scalar v_d2_div = forall<double>(v_d2_div_calc,v_d2_ahead,v_d2_behind,v_d2_behind_dp,v_d2_behind_dm);
					h_lambda += v_d2_div;
				}
			}
			Scalar Fast_MS_speed = forall<double>(Fast_MS_speed_calc, W_bar, d, gamma);
			Scalar Fast_MS_speed_behind = alias(Fast_MS_speed,Point::Basis(d)*(1));
			Scalar Fast_MS_speed_min = forall<double>(Fast_MS_speed_min_calc,Fast_MS_speed,Fast_MS_speed_behind);
			Scalar Visc_coef = forall<double>(Visc_coef_calc,h_lambda,Fast_MS_speed_min);
			Vector a_U_behind = m_behind_shift[d](a_U);
			Vector mu_f = forall<double,NUMCOMPS>(mu_f_calc, Visc_coef, a_U, a_U_behind);
			a_Rhs += m_divergence[d](mu_f);
	    }
		a_Rhs *= -1./a_dx;
	}
  
  
  
}
