#include "Proto.H"
#include "MHD_Mapping.H"
#include "CommonTemplates.H"
#include "Proto_Timer.H"
#include "Proto_WriteBoxData.H"

const double C1_fix = 0.00;
const double r_in = 0.3;
const double r_out = 0.6;
extern int grid_type_global;
typedef BoxData<double,1,1,1> Scalar;
typedef BoxData<double,NUMCOMPS,1,1> Vector;

namespace MHD_Mapping{
	
    PROTO_KERNEL_START
	void eta_to_xF(Var<double,DIM>& a_x,
						const Var<double,DIM>& a_eta)
	{
		
#if DIM == 2		
		if (grid_type_global == 0){
			a_x(0) = a_eta(0);
			a_x(1) = a_eta(1);
		}
		if (grid_type_global == 1){
			double C1 = C1_fix;
			double C2 = C1;
			a_x(0) = a_eta(0) + C1*sin(2.0*PI*a_eta(0))*sin(2.0*PI*a_eta(1));
			a_x(1) = a_eta(1) + C2*sin(2.0*PI*a_eta(0))*sin(2.0*PI*a_eta(1));
		}
		if (grid_type_global == 2){
			a_x(0) = (r_in + a_eta(0)*(r_out-r_in))*cos(2.0*PI*a_eta(1));
			a_x(1) = (r_in + a_eta(0)*(r_out-r_in))*sin(2.0*PI*a_eta(1));
		}
#endif		
	}
	PROTO_KERNEL_END(eta_to_xF, eta_to_x)
	
	
	void eta_to_x_calc(BoxData<double,DIM>& a_x,
					   const BoxData<double,DIM>& a_eta)
	{
		a_x = forall<double,DIM>(eta_to_x,a_eta);
	}
	
	
	
	PROTO_KERNEL_START
	void dot_pro_calcXF(Var<double,1>& a_dot_pro,
					 const Var<double,1>& a_d_perp_N_s,
					 const Var<double,1>& a_d_perp_X)
	{
			a_dot_pro(0) = (a_d_perp_N_s(0)*a_d_perp_X(0));			
	}
	PROTO_KERNEL_END(dot_pro_calcXF, dot_pro_calcX)
	
	
	PROTO_KERNEL_START
	void X_f_mapped1D_calcF(Var<double,1>& a_X_f_mapped1D,
					 const Var<double,1>& a_X_ave_f,
					 const Var<double,1>& a_N_s_d_ave_f,
					 const Var<double,1>& a_dot_pro_sum)
	{
			a_X_f_mapped1D(0) = a_N_s_d_ave_f(0)*a_X_ave_f(0) + a_dot_pro_sum(0)/12.0;			
	}
	PROTO_KERNEL_END(X_f_mapped1D_calcF, X_f_mapped1D_calc)

	
	
	
	PROTO_KERNEL_START
	void X_ave_f_calcF( const Point& a_pt,
					 Var<double,1>& a_X_ave_f,
					 int a_s,
					 int a_d,
					 const double a_dx)
	{	
		double x_loc = a_pt[0];
		double y_loc = a_pt[1];	
		double x_hi = (a_pt[0] + 1.0)*a_dx;
		double x_lo = (a_pt[0] - 0.0)*a_dx;
		double y_hi = (a_pt[1] + 1.0)*a_dx;
		double y_lo = (a_pt[1] - 0.0)*a_dx;
#if DIM == 2		
		if (grid_type_global == 0){ 
			double C1 = 0.0;
			double C2 = C1;
			if (a_s == 0 && a_d == 0) a_X_ave_f(0) = ((x_lo*y_hi)
													 -(x_lo*y_lo))/DIM/a_dx;
			if (a_s == 1 && a_d == 0) a_X_ave_f(0) = ((y_hi*y_hi/2.0)
													 -(y_lo*y_lo/2.0))/DIM/a_dx;
			if (a_s == 0 && a_d == 1) a_X_ave_f(0) = ((x_hi*x_hi/2.0)
													 -(x_lo*x_lo/2.0))/DIM/a_dx;
			if (a_s == 1 && a_d == 1) a_X_ave_f(0) = ((x_hi*y_lo)
													 -(x_lo*y_lo))/DIM/a_dx;
		}
		
		if (grid_type_global == 1){ 
			double C1 = C1_fix;
			double C2 = C1;
			if (a_s == 0 && a_d == 0) a_X_ave_f(0) = ((x_lo*y_hi-(C1/2./PI)*sin(2.*PI*x_lo)*cos(2.*PI*y_hi))
													 -(x_lo*y_lo-(C1/2./PI)*sin(2.*PI*x_lo)*cos(2.*PI*y_lo)))/DIM/a_dx;
			if (a_s == 1 && a_d == 0) a_X_ave_f(0) = ((y_hi*y_hi/2.0-(C2/2./PI)*sin(2.*PI*x_lo)*cos(2.*PI*y_hi))
													 -(y_lo*y_lo/2.0-(C2/2./PI)*sin(2.*PI*x_lo)*cos(2.*PI*y_lo)))/DIM/a_dx;
			if (a_s == 0 && a_d == 1) a_X_ave_f(0) = ((x_hi*x_hi/2.0-(C1/2./PI)*cos(2.*PI*x_hi)*sin(2.*PI*y_lo))
													 -(x_lo*x_lo/2.0-(C1/2./PI)*cos(2.*PI*x_lo)*sin(2.*PI*y_lo)))/DIM/a_dx;
			if (a_s == 1 && a_d == 1) a_X_ave_f(0) = ((x_hi*y_lo-(C2/2./PI)*cos(2.*PI*x_hi)*sin(2.*PI*y_lo))
													 -(x_lo*y_lo-(C2/2./PI)*cos(2.*PI*x_lo)*sin(2.*PI*y_lo)))/DIM/a_dx;
		}
		
		if (grid_type_global == 2){ 
			double r_diff = r_out-r_in;
			if (a_s == 0 && a_d == 0) a_X_ave_f(0) = (((r_in+x_lo*r_diff)*sin(2.*PI*y_hi)/2./PI)
													 -((r_in+x_lo*r_diff)*sin(2.*PI*y_lo)/2./PI))/DIM/a_dx;
			if (a_s == 1 && a_d == 0) a_X_ave_f(0) = -(((r_in+x_lo*r_diff)*cos(2.*PI*y_hi)/2./PI)
													  -((r_in+x_lo*r_diff)*cos(2.*PI*y_lo)/2./PI))/DIM/a_dx;
			if (a_s == 0 && a_d == 1) a_X_ave_f(0) = (((r_in*x_hi+x_hi*x_hi*r_diff/2.0)*cos(2.*PI*y_lo))
													 -((r_in*x_lo+x_lo*x_lo*r_diff/2.0)*cos(2.*PI*y_lo)))/DIM/a_dx;
			if (a_s == 1 && a_d == 1) a_X_ave_f(0) = (((r_in*x_hi+x_hi*x_hi*r_diff/2.0)*sin(2.*PI*y_lo))
													 -((r_in*x_lo+x_lo*x_lo*r_diff/2.0)*sin(2.*PI*y_lo)))/DIM/a_dx;
		}
#endif
		

	  
	}
	PROTO_KERNEL_END(X_ave_f_calcF, X_ave_f_calc)
	
	
	
	PROTO_KERNEL_START
	void N_ave_f_calcF( const Point& a_pt,
					 Var<double,1>& a_N_ave_f,
					 int a_s,
					 int a_d,
					 const double a_dx)
	{	
		double x_loc = a_pt[0];
		double y_loc = a_pt[1];	
		double x_hi = (a_pt[0] + 1.0)*a_dx;
		double x_lo = (a_pt[0] - 0.0)*a_dx;
		double y_hi = (a_pt[1] + 1.0)*a_dx;
		double y_lo = (a_pt[1] - 0.0)*a_dx;
		
		if (grid_type_global == 0){
			double C1 = 0.0;
			double C2 = C1;
			if (a_s == 0 && a_d == 0) a_N_ave_f(0) = 1.0;
			if (a_s == 1 && a_d == 0) a_N_ave_f(0) = 0.0;
			if (a_s == 0 && a_d == 1) a_N_ave_f(0) = 0.0;
			if (a_s == 1 && a_d == 1) a_N_ave_f(0) = 1.0;	
		}
		
		if (grid_type_global == 1){
			double C1 = C1_fix;
			double C2 = C1;
			if (a_s == 0 && a_d == 0) a_N_ave_f(0) = ((y_hi+C2*sin(2.*PI*x_lo)*sin(2.*PI*y_hi))
													 -(y_lo+C2*sin(2.*PI*x_lo)*sin(2.*PI*y_lo)))/a_dx;
			if (a_s == 1 && a_d == 0) a_N_ave_f(0) = -((x_lo+C1*sin(2.*PI*x_lo)*sin(2.*PI*y_hi))
													  -(x_lo+C1*sin(2.*PI*x_lo)*sin(2.*PI*y_lo)))/a_dx;
			if (a_s == 0 && a_d == 1) a_N_ave_f(0) = -((y_lo+C2*sin(2.*PI*x_hi)*sin(2.*PI*y_lo))
													  -(y_lo+C2*sin(2.*PI*x_lo)*sin(2.*PI*y_lo)))/a_dx;
			if (a_s == 1 && a_d == 1) a_N_ave_f(0) = ((x_hi+C1*sin(2.*PI*x_hi)*sin(2.*PI*y_lo))
													 -(x_lo+C1*sin(2.*PI*x_lo)*sin(2.*PI*y_lo)))/a_dx;		
		}

		if (grid_type_global == 2){
			double r_diff = r_out-r_in;
			if (a_s == 0 && a_d == 0) a_N_ave_f(0) = (((r_in + x_lo*r_diff)*sin(2.*PI*y_hi))
													 -((r_in + x_lo*r_diff)*sin(2.*PI*y_lo)))/a_dx;
			if (a_s == 1 && a_d == 0) a_N_ave_f(0) = -(((r_in + x_lo*r_diff)*cos(2.*PI*y_hi))
													  -((r_in + x_lo*r_diff)*cos(2.*PI*y_lo)))/a_dx;
			if (a_s == 0 && a_d == 1) a_N_ave_f(0) = -(((r_in + x_hi*r_diff)*sin(2.*PI*y_lo))
													  -((r_in + x_lo*r_diff)*sin(2.*PI*y_lo)))/a_dx;
			if (a_s == 1 && a_d == 1) a_N_ave_f(0) = (((r_in + x_hi*r_diff)*cos(2.*PI*y_lo))
													 -((r_in + x_lo*r_diff)*cos(2.*PI*y_lo)))/a_dx;		
		}

		
	}
	PROTO_KERNEL_END(N_ave_f_calcF, N_ave_f_calc)
	
	
	PROTO_KERNEL_START
	void d_U_calcF(  State& a_d_U,
					 const State& a_UJ_ahead,
					 const State& a_UJ_behind,
					 const Var<double,1>& a_J_ahead,
					 const Var<double,1>& a_J_behind)
	{	
		for (int i=0; i< NUMCOMPS; i++){
			a_d_U(i) = ((a_UJ_ahead(i)/a_J_ahead(0))-(a_UJ_behind(i)/a_J_behind(0)))/(2.0);			
		}

	}
	PROTO_KERNEL_END(d_U_calcF, d_U_calc)
	
	
	
	PROTO_KERNEL_START
	void a_U_demapped_calcF(  State& a_U_demapped,
					 const Var<double,1>& a_Jacobian_ave,
					 const State& a_a_U,
					 const State& a_dot_pro_sum2)
	{	
		for (int i=0; i< NUMCOMPS; i++){
			a_U_demapped(i) = (a_a_U(i)-(a_dot_pro_sum2(i)/12.0))/a_Jacobian_ave(0);			
		}

	}
	PROTO_KERNEL_END(a_U_demapped_calcF, a_U_demapped_calc)
	
	
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
	
	
	void N_ave_f_calc_func(BoxData<double,1>& a_N_ave_f,
						   const int a_s,
						   const int a_d,
						   const double a_dx)
	{
		forallInPlace_p(N_ave_f_calc, a_N_ave_f, a_s, a_d, a_dx);
	}							   
	
	
	void Jacobian_Ave_calc(BoxData<double,1>& a_Jacobian_ave,
						   const double a_dx,
						   const Box& a_dbx0)
	{
		
		
		static Stencil<double> m_divergence[DIM];
		static Stencil<double> m_derivative[DIM];
		static bool initialized = false;
		if(!initialized)
		{
		  for (int dir = 0; dir < DIM; dir++)
		  {
			m_divergence[dir] = Stencil<double>::FluxDivergence(dir);
			m_derivative[dir] = Stencil<double>::Derivative(1,dir,2);
		  }
		  initialized =  true;
		}
		

		Scalar N_s_d_ave_f(a_dbx0), X_ave_f(a_dbx0);

		a_Jacobian_ave.setVal(0.0);
		
		for (int d = 0; d < DIM; d++){
			Scalar X_f_mapped(a_dbx0);
			X_f_mapped.setVal(0.0);
			for (int s = 0; s < DIM; s++){
				forallInPlace_p(N_ave_f_calc, N_s_d_ave_f, s, d, a_dx);
				forallInPlace_p(X_ave_f_calc, X_ave_f, s, d, a_dx);
				Scalar dot_pro_sum(a_dbx0);
				dot_pro_sum.setVal(0.0);
				for (int s_temp = 0; s_temp < DIM; s_temp++){
					if (s_temp != d){
						Scalar d_perp_N_s = m_derivative[s_temp](N_s_d_ave_f);
						Scalar d_perp_X = m_derivative[s_temp](X_ave_f);
						Scalar dot_pro = forall<double,1>(dot_pro_calcX,d_perp_N_s,d_perp_X);
						dot_pro_sum += dot_pro;
					}
				}
				Scalar X_f_mapped1D = forall<double,1>(X_f_mapped1D_calc,X_ave_f,N_s_d_ave_f,dot_pro_sum);
				X_f_mapped += X_f_mapped1D;
			}
			a_Jacobian_ave += m_divergence[d](X_f_mapped);
		}
		a_Jacobian_ave *= 1./a_dx;
		
	}
	
	
	
	void JU_to_U_calc(BoxData<double,NUMCOMPS>& a_U_demapped,
				      const BoxData<double,NUMCOMPS>& a_a_U,
					  const BoxData<double,1>& a_Jacobian_ave,
				      const Box& a_dbx0)
	{
		

		static Stencil<double> m_ahead_shift[DIM];
		static Stencil<double> m_behind_shift[DIM];
		static Stencil<double> m_derivative[DIM];
		static bool initialized = false;
		if(!initialized)
		{
		  for (int dir = 0; dir < DIM; dir++)
		  {
			m_ahead_shift[dir] = 1.0*Shift(Point::Basis(dir)*(1));
			m_behind_shift[dir] = 1.0*Shift(Point::Basis(dir)*(-1));
			m_derivative[dir] = Stencil<double>::Derivative(1,dir,2);
		  }
		  initialized =  true;
		}
		
		
		
		Vector dot_pro_sum2(a_dbx0);
		dot_pro_sum2.setVal(0.0);
		for (int d = 0; d < DIM; d++){
			Vector UJ_ahead = m_ahead_shift[d](a_a_U);
			Vector UJ_behind = m_behind_shift[d](a_a_U);
			Scalar J_ahead = m_ahead_shift[d](a_Jacobian_ave);
			Scalar J_behind = m_behind_shift[d](a_Jacobian_ave);
			Vector d_U = forall<double,NUMCOMPS>(d_U_calc,UJ_ahead,UJ_behind,J_ahead,J_behind);
			Scalar d_J = m_derivative[d](a_Jacobian_ave);
			Vector dot_pro = forall<double,NUMCOMPS>(dot_pro_calcF,d_J,d_U);
			dot_pro_sum2 += dot_pro;
		}

		a_U_demapped = forall<double,NUMCOMPS>(a_U_demapped_calc,a_Jacobian_ave,a_a_U,dot_pro_sum2);
	}
						  
	
	
	
	
							   
}