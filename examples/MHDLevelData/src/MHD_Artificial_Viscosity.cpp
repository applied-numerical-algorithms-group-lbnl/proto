#include "Proto.H"
#include "MHDOp.H"
#include "MHD_Mapping.H"
#include "MHD_Artificial_Viscosity.H"

typedef BoxData<double,1,1,1> Scalar;
typedef BoxData<double,NUMCOMPS,1,1> Vector;

namespace MHD_Artificial_Viscosity {
	
	
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
	
	
	
	

// Used to implement artificial viscosity
	void step(BoxData<double,NUMCOMPS>& a_Rhs,
            const BoxData<double,NUMCOMPS>& a_U,
            const Box& a_rangeBox,
            const double a_dx,
            const double a_gamma,
            Reduction<double>& a_Rxn,
            bool a_computeMaxWaveSpeed,
            bool a_callBCs)
	{

	    Box dbx0 = a_U.box();
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
		
		Scalar Jacobian_ave(dbx0);
		MHD_Mapping::Jacobian_Ave_calc(Jacobian_ave,a_dx,dbx0);
		
		Vector a_U_demapped(dbx0);
		MHD_Mapping::JU_to_U_calc(a_U_demapped, a_U, Jacobian_ave, dbx0);
		
		Vector W_bar(dbx0);
		MHDOp::consToPrimcalc(W_bar,a_U_demapped,gamma);
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
			Vector a_U_behind = m_behind_shift[d](a_U_demapped);
			Vector mu_f = forall<double,NUMCOMPS>(mu_f_calc, Visc_coef, a_U_demapped, a_U_behind);
			a_Rhs += m_divergence[d](mu_f);
	    }
		a_Rhs *= -1./a_dx;
	}
}