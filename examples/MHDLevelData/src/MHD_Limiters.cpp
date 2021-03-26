#include "Proto.H"
#include "MHD_Limiters.H"
#include "CommonTemplates.H"
#include "Proto_Timer.H"
#include "Proto_WriteBoxData.H"

typedef BoxData<double,1,1,1> Scalar;
typedef BoxData<double,NUMCOMPS,1,1> Vector;

namespace MHD_Limiters{
	
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
	void W_ave_diff_calcF(State& W_ave_diff,
					 const State& a_W_ave,
					 const State& a_W_ahead)
	{
		for (int i=0; i< NUMCOMPS; i++){
			W_ave_diff(i) = a_W_ave(i) - a_W_ahead(i);
		}
	}
	PROTO_KERNEL_END(W_ave_diff_calcF, W_ave_diff_calc)

	
	

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
	void del3_W_C_calcF(State& a_del3_W_C,
					 const State& a_del2_W_c_ahead,
					 const State& a_del2_W_c_behind)
	{
		for (int i=0; i< NUMCOMPS; i++){
			a_del3_W_C(i) = (a_del2_W_c_ahead(i) - a_del2_W_c_behind(i))/2.0;
		}
	}
	PROTO_KERNEL_END(del3_W_C_calcF, del3_W_C_calc)

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

    const bool printlim = false;

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
					 const State& a_del3_W_L,
					 const State& a_del3_W_R,
					 const State& a_del3_W_C,
					 int   a_dir,
					 const double a_dx)
	{
		double rhs = 1.0 - 1.0e-12;
		double rhs_test, lhs_test, rhs_test2, lhs_test2,rhs_test_a, lhs_test_a;
		double x_loc = a_pt[0];//*a_dx + 0.5*a_dx;
		double y_loc = a_pt[1];//*a_dx + 0.5*a_dx;
		//cerr << x_loc << " " << y_loc << endl;
		double x_loc2 = a_pt[0]*a_dx + 0.5*a_dx;
		double y_loc2 = a_pt[1]*a_dx + 0.5*a_dx;
		
		for (int i=0; i< NUMCOMPS; i++)
		{ 
	      if (i != (NUMCOMPS + a_dir - 2))
	      //if (i != 100)
		  {	  
	        if ((a_del_W_f_m(i) * a_del_W_f_p(i) <= 0.0) || ((a_W_ave(i) - a_W_ave_behind2(i))*(a_W_ave_ahead2(i) - a_W_ave(i)) <= 0.0)){
				if (a_rho_i(i) < rhs){
					rhs_test = a_del3_W_max(i) - a_del3_W_min(i);
					lhs_test = 0.1*std::max({std::abs(a_del3_W_min(i)), std::abs(a_del3_W_max(i))});
					double prevval = a_W_ave_low_ahead(i);
					lhs_test_a = std::min({std::abs(a_del3_W_L(i)),std::abs(a_del3_W_R(i))})*1.1;
					rhs_test_a = std::abs(a_del3_W_C(i));
					if (lhs_test <= rhs_test){
					//if (!(lhs_test_a > rhs_test_a && a_del3_W_L(i)*a_del3_W_R(i) > 0.)){
						rhs_test2 = 2.0*std::abs(a_del_W_f_m(i));
						lhs_test2 = std::abs(a_del_W_f_p(i));
						if (a_del_W_f_m(i) * a_del_W_f_p(i) < 0.0){
							a_W_ave_low_ahead_limited(i) = a_W_ave(i) + a_rho_i(i)*a_del_W_f_p(i);
							if (printlim){
								if (std::abs(a_W_ave_low_ahead_limited(i) - a_W_ave_low_ahead(i)) > 1.0e-12 && x_loc2 > 0.0 && x_loc2 < 1.0 && y_loc2 > 0.0 && y_loc2 < 1.0) cout << "x = " << x_loc << " y = " << y_loc << " Component = "<< i << " dir = " << a_dir <<" Condition : 1_a_L) Before limit: " << a_W_ave_low_ahead(i) << " Change: " << a_W_ave_low_ahead_limited(i) - a_W_ave_low_ahead(i) << " rho: " << a_rho_i(i) << endl;
								//if (std::abs(a_W_ave_low_ahead_limited(i) - a_W_ave_low_ahead(i)) > 1.0e-14 && x_loc2 > 0.46 && x_loc2 < .54 && y_loc2 > 0.35 && y_loc2 < .41) cout << lhs_test << " " << rhs_test << endl;
							}
						} else if (lhs_test2 >= rhs_test2){
							a_W_ave_low_ahead_limited(i) = a_W_ave(i) + 2.0*(1.0-a_rho_i(i))*a_del_W_f_m(i) + a_rho_i(i)*a_del_W_f_p(i);
							if (printlim){
								if (std::abs(a_W_ave_low_ahead_limited(i) - a_W_ave_low_ahead(i)) > 1.0e-12 && x_loc2 > 0.0 && x_loc2 < 1.0 && y_loc2 > 0.0 && y_loc2 < 1.0) cout << "x = " << x_loc << " y = " << y_loc << " Component = "<< i << " dir = " << a_dir <<" Condition : 1_c_L) Before limit: " << a_W_ave_low_ahead(i) << " Change: " << a_W_ave_low_ahead_limited(i) - a_W_ave_low_ahead(i) <<" rho: " << a_rho_i(i)<< endl;
								//if (std::abs(a_W_ave_low_ahead_limited(i) - a_W_ave_low_ahead(i)) > 1.0e-14 && x_loc2 > 0.46 && x_loc2 < .54 && y_loc2 > 0.35 && y_loc2 < .41) cout << lhs_test << " " << rhs_test << endl;
							}
						}
					}
					//}
				}
			} else {
				rhs_test2 = 2.0*std::abs(a_del_W_f_m(i));
				lhs_test2 = std::abs(a_del_W_f_p(i));
				if (lhs_test2 >= rhs_test2){
					a_W_ave_low_ahead_limited(i) = a_W_ave(i) + 2.0*a_del_W_f_m(i);
					if (printlim){
						if (std::abs(a_W_ave_low_ahead_limited(i) - a_W_ave_low_ahead(i)) > 1.0e-12 && x_loc2 > 0.0 && x_loc2 < 1.0 && y_loc2 > 0.0 && y_loc2 < 1.0) cout << "x = " << x_loc << " y = " << y_loc << " Component = "<< i << " dir = " << a_dir <<" Condition : 2_b_L) Before limit: " << a_W_ave_low_ahead(i) << " Change: " << a_W_ave_low_ahead_limited(i) - a_W_ave_low_ahead(i) << " rho: " << a_rho_i(i) << endl;
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
					 const State& a_del3_W_L,
					 const State& a_del3_W_R,
					 const State& a_del3_W_C,
					 int   a_dir,
					 const double a_dx)
	{
		double rhs = 1.0 - 1.0e-12;
		double rhs_test, lhs_test, rhs_test2, lhs_test2,rhs_test_a, lhs_test_a;
		double x_loc = a_pt[0];//*a_dx + 0.5*a_dx;
		double y_loc = a_pt[1];//*a_dx + 0.5*a_dx;
		double x_loc2 = a_pt[0]*a_dx + 0.5*a_dx;
		double y_loc2 = a_pt[1]*a_dx + 0.5*a_dx;
		
		
		for (int i=0; i< NUMCOMPS; i++)
		{
		  if (i != (NUMCOMPS + a_dir - 2))
	      //if (i != 100)
		  {	
	        if ((a_del_W_f_m(i) * a_del_W_f_p(i) <= 0.0) || ((a_W_ave(i) - a_W_ave_behind2(i))*(a_W_ave_ahead2(i) - a_W_ave(i)) <= 0.0)){
				if (a_rho_i(i) < rhs){
					rhs_test = a_del3_W_max(i) - a_del3_W_min(i);
					lhs_test = 0.1*std::max({std::abs(a_del3_W_min(i)), std::abs(a_del3_W_max(i))});
					lhs_test_a = std::min({std::abs(a_del3_W_L(i)),std::abs(a_del3_W_R(i))})*1.1;
					rhs_test_a = std::abs(a_del3_W_C(i));
					if (lhs_test <= rhs_test){
					//if (!(lhs_test_a > rhs_test_a && a_del3_W_L(i)*a_del3_W_R(i) > 0.)){
						
						rhs_test2 = 2.0*std::abs(a_del_W_f_p(i));
						lhs_test2 = std::abs(a_del_W_f_m(i));
						if (a_del_W_f_m(i) * a_del_W_f_p(i) < 0.0){
							a_W_ave_high_limited(i) = a_W_ave(i) - a_rho_i(i)*a_del_W_f_m(i);
							if (printlim){
								if (std::abs(a_W_ave_high_limited(i) - a_W_ave_high(i)) > 1.0e-12 && x_loc2 > 0.0 && x_loc2 < 1.0 && y_loc2 > 0.0 && y_loc2 < 1.0) cout << "x = " << x_loc << " y = " << y_loc << " Component = "<< i << " dir = " << a_dir <<" Condition : 1_a_R) Before limit: " << a_W_ave_high(i) << " Change: " << a_W_ave_high_limited(i) - a_W_ave_high(i) << " rho: " << a_rho_i(i) << endl;
								//if (std::abs(a_W_ave_high_limited(i) - a_W_ave_high(i)) > 1.0e-14 && x_loc2 > 0.46 && x_loc2 < .54 && y_loc2 > 0.35 && y_loc2 < .41) cout << lhs_test << " " << rhs_test << endl;
							}
						} else if (lhs_test2 >= rhs_test2){
							a_W_ave_high_limited(i) = a_W_ave(i) - 2.0*(1.0-a_rho_i(i))*a_del_W_f_p(i) - a_rho_i(i)*a_del_W_f_m(i);
							if (printlim){
								if (std::abs(a_W_ave_high_limited(i) - a_W_ave_high(i)) > 1.0e-12 && x_loc2 > 0.0 && x_loc2 < 1.0 && y_loc2 > 0.0 && y_loc2 < 1.0) cout << "x = " << x_loc << " y = " << y_loc << " Component = "<< i << " dir = " << a_dir <<" Condition : 1_b_R) Before limit: " << a_W_ave_high(i) << " Change: " << a_W_ave_high_limited(i) - a_W_ave_high(i) << " rho: " << a_rho_i(i) << " " << lhs_test2 << " " << rhs_test2/2.0 << endl;
								//if (std::abs(a_W_ave_high_limited(i) - a_W_ave_high(i)) > 1.0e-14 && x_loc2 > 0.46 && x_loc2 < .54 && y_loc2 > 0.35 && y_loc2 < .41) cout << lhs_test << " " << rhs_test << endl;
							}
						}
					}
					//}
				}
			} else {
				rhs_test2 = 2.0*std::abs(a_del_W_f_p(i));
				lhs_test2 = std::abs(a_del_W_f_m(i));
				if (lhs_test2 >= rhs_test2){
					a_W_ave_high_limited(i) = a_W_ave(i) - 2.0*a_del_W_f_p(i);
					if (printlim){
						if (std::abs(a_W_ave_high_limited(i) - a_W_ave_high(i)) > 1.0e-12 && x_loc2 > 0.0 && x_loc2 < 1.0 && y_loc2 > 0.0 && y_loc2 < 1.0) cout << "x = " << x_loc << " y = " << y_loc << " Component = "<< i << " dir = " << a_dir <<" Condition : 2_a_R) Before limit: " << a_W_ave_high(i) << " Change: " << a_W_ave_high_limited(i) - a_W_ave_high(i) << " rho: " << a_rho_i(i) << endl;
					}
				}
			}		
		  }
		}
	}
	PROTO_KERNEL_END(limiter_high_calcF, limiter_high_calc)


	
	
	
	
	PROTO_KERNEL_START
	void print_limiter_data_calcF(
	                 const Point& a_pt, 
	                 State& W_ave_high_temp,
					 const State& a_W_ave,
					 const State& a_W_ave_low_ahead,
					 const State& a_W_ave_high,
					 const State& a_del_W_f_p,
					 const State& a_del_W_f_m,
					 const State& a_del2_W_c,
					 const State& a_del2_W_f,
					 const State& a_del3_W,
					 const State& a_del2_W_lim,
					 const State& a_rho_i,
					 const State& a_del3_W_min,
					 const State& a_del3_W_max,
					 const State& a_del3_W_L,
					 const State& a_del3_W_R,
					 const State& a_del3_W_C,
					 int   a_dir,
					 const double a_dx)
	{
		
		double x_loc = a_pt[0];//*a_dx + 0.5*a_dx;
		double y_loc = a_pt[1];//*a_dx + 0.5*a_dx;
		double x_loc2 = a_pt[0]*a_dx + 0.5*a_dx;
		double y_loc2 = a_pt[1]*a_dx + 0.5*a_dx;
		
		
		//for (int i=0; i< NUMCOMPS; i++)
		for (int i=3; i<= 3; i++)
		{
			if (printlim){
				if (x_loc >= 1 && x_loc <= 10 && y_loc == 9 ) {
					cout << setw(20) << setprecision(10) << x_loc 
					     << setw(20) << setprecision(10) << y_loc 
						 << setw(20) << setprecision(10) << a_W_ave(i) 
					     << setw(20) << setprecision(10) << a_W_ave_low_ahead(i) 
					     << setw(20) << setprecision(10) << a_W_ave_high(i) 
					     << setw(20) << setprecision(10) << a_del_W_f_p(i) 
					     << setw(20) << setprecision(10) << a_del_W_f_m(i) 
					     << setw(20) << setprecision(10) << a_del2_W_c(i) 
					     << setw(20) << setprecision(10) << a_del2_W_f(i) 
					     << setw(20) << setprecision(10) << a_del3_W(i) 
					     << setw(20) << setprecision(10) << a_del2_W_lim(i) 
					     << setw(20) << setprecision(10) << a_rho_i(i) 
					     << setw(20) << setprecision(10) << a_del3_W_min(i)
					     << setw(20) << setprecision(10) << a_del3_W_max(i) 
					     << setw(20) << setprecision(10) << a_del3_W_L(i) 
					     << setw(20) << setprecision(10) << a_del3_W_R(i) 
					     << setw(20) << setprecision(10) << a_del3_W_C(i) 
						 << endl;
				}
			}
								
		  
		}
	}
	PROTO_KERNEL_END(print_limiter_data_calcF, print_limiter_data_calc)
	
	
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
	
	
	void MHD_Limiters(BoxData<double,NUMCOMPS>& a_W_ave_low_lim_flat,
				  BoxData<double,NUMCOMPS>& a_W_ave_high_lim_flat,
				  const BoxData<double,NUMCOMPS>& a_W_ave_low,
				  const BoxData<double,NUMCOMPS>& a_W_ave_high,
				  const BoxData<double,NUMCOMPS>& a_W_ave,
				  const int a_d,
                  const double a_dx)
	{
		bool out_debug_data = false;
		static Stencil<double> m_copy;
		static Stencil<double> m_ahead_shift[DIM];
		static Stencil<double> m_behind_shift[DIM];
		static Stencil<double> m_copy_f[DIM];
		static bool initialized = false;
		if(!initialized)
		{
		  m_copy = 1.0*Shift(Point::Zeros());
		  for (int dir = 0; dir < DIM; dir++)
		  {
			m_ahead_shift[dir] = 1.0*Shift(Point::Basis(dir)*(1));
			m_behind_shift[dir] = 1.0*Shift(Point::Basis(dir)*(-1));
			m_copy_f[dir] = 1.0*Shift(Point::Zeros());
		  }
		  initialized =  true;
		}
		//Limiter Starts here
			//Vector W_ave_low_ahead = alias(a_W_ave_low,Point::Basis(a_d)*(-1));
			Vector W_ave_low_ahead = m_ahead_shift[a_d](a_W_ave_low);
			Vector del_W_f_m = forall<double,NUMCOMPS>(del_W_f_m_calc, a_W_ave, a_W_ave_high);
			Vector del_W_f_p = forall<double,NUMCOMPS>(del_W_f_p_calc, a_W_ave, W_ave_low_ahead);
			Vector del2_W_f  = forall<double,NUMCOMPS>(del2_W_f_calc, a_W_ave, a_W_ave_high, W_ave_low_ahead);
			//Vector W_ave_ahead = alias(a_W_ave,Point::Basis(a_d)*(-1));
			Vector W_ave_ahead = m_ahead_shift[a_d](a_W_ave);
			Vector W_ave_diff = forall<double,NUMCOMPS>(W_ave_diff_calc, a_W_ave, W_ave_ahead);
			//Vector W_ave_ahead2 = alias(a_W_ave,Point::Basis(a_d)*(-2));
			Vector W_ave_ahead2 = m_ahead_shift[a_d](W_ave_ahead);
			//Vector W_ave_behind = alias(a_W_ave,Point::Basis(a_d)*(1));
			Vector W_ave_behind = m_behind_shift[a_d](a_W_ave);
			//Vector W_ave_behind2 = alias(a_W_ave,Point::Basis(a_d)*(2));
			Vector W_ave_behind2 = m_behind_shift[a_d](W_ave_behind);
			Vector del2_W_c  = forall<double,NUMCOMPS>( del2_W_c_calc, a_W_ave, W_ave_behind, W_ave_ahead);
			//Vector del2_W_c_ahead = alias(del2_W_c,Point::Basis(a_d)*(-1));
			Vector del2_W_c_ahead = m_ahead_shift[a_d](del2_W_c);
			//Vector del2_W_c_behind = alias(del2_W_c,Point::Basis(a_d)*(1));
			Vector del2_W_c_behind = m_behind_shift[a_d](del2_W_c);
			Vector del3_W = forall<double,NUMCOMPS>( del3_W_calc, del2_W_c, del2_W_c_behind);
			Vector del3_W_C = forall<double,NUMCOMPS>( del3_W_C_calc, del2_W_c_ahead, del2_W_c_behind);
			Vector del3_W_L = forall<double,NUMCOMPS>( del3_W_calc, del2_W_c, del2_W_c_behind);
			Vector del3_W_R = forall<double,NUMCOMPS>( del3_W_calc, del2_W_c_ahead, del2_W_c);
			//Vector del3_W_behind = alias(del3_W,Point::Basis(a_d)*(1));
			Vector del3_W_behind = m_behind_shift[a_d](del3_W);
			//Vector del3_W_ahead = alias(del3_W,Point::Basis(a_d)*(-1));
			Vector del3_W_ahead = m_ahead_shift[a_d](del3_W);
			//Vector del3_W_ahead2 = alias(del3_W,Point::Basis(a_d)*(-2));
			Vector del3_W_ahead2 = m_ahead_shift[a_d](del3_W_ahead);
			Vector del2_W_lim = forall<double,NUMCOMPS>(del2_W_lim_calc, del2_W_f, del2_W_c, del2_W_c_ahead, del2_W_c_behind);
			Vector rho_i = forall<double,NUMCOMPS>( rho_i_calc, del2_W_f, del2_W_lim, a_W_ave, 
													W_ave_ahead, W_ave_ahead2, W_ave_behind, W_ave_behind2);
			Vector del3_W_min = forall<double,NUMCOMPS>(del3_W_min_calc, del3_W_behind, 
														 del3_W, del3_W_ahead, del3_W_ahead2);	
			Vector del3_W_max = forall<double,NUMCOMPS>( del3_W_max_calc, del3_W_behind, 
														 del3_W, del3_W_ahead, del3_W_ahead2);				 
			Vector W_ave_low_ahead_limited = m_copy_f[a_d](W_ave_low_ahead);
			Vector W_ave_high_limited = m_copy_f[a_d](a_W_ave_high);
            forallInPlace_p( limiter_low_calc,W_ave_low_ahead_limited, W_ave_low_ahead, del_W_f_m, del_W_f_p, a_W_ave, W_ave_ahead2, W_ave_behind2, rho_i, del3_W_max, del3_W_min,del3_W_L,del3_W_R,del3_W_C, a_d, a_dx);
			forallInPlace_p( limiter_high_calc, W_ave_high_limited, a_W_ave_high, del_W_f_m, del_W_f_p, a_W_ave, W_ave_ahead2, W_ave_behind2, rho_i, del3_W_max, del3_W_min,del3_W_L,del3_W_R,del3_W_C, a_d, a_dx);
			
			
			// if (a_d==0 && out_debug_data){
				// std::string filename="W_ave";
			    // // const char* vectNames[6];
				// // vectNames[0] = "rho";
				// // vectNames[1] = "momentum_x";
				// // vectNames[2] = "momentum_y";
				// // vectNames[3] = "energy";
				// // vectNames[4] = "B_x";
				// // vectNames[5] = "B_y";
				// // double origin[DIM];
				// // for (int ii = 0; ii < DIM; ii++)
				// // {
						// // origin[ii] = 0.0;
				// // }
				// // WriteBoxData(filename.c_str(),W_ave,vectNames,origin,a_dx);
				
				// //copyTo(W_ave,a_rangeBox);
				// // Vector outBD = m_copy(W_ave,a_rangeBox);
				// // BoxData<double> W_ave_out = slice(outBD,1);
				// // WriteBoxData(filename.c_str(),W_ave_out,a_dx);
				
				// Vector W_ave_high_temp = m_copy_f[d](W_ave_high);
				// forallInPlace_p( print_limiter_data_calc,W_ave_high_temp,W_ave,W_ave_low_ahead,W_ave_high,del_W_f_p,del_W_f_m,del2_W_c,del2_W_f,del3_W,del2_W_lim,rho_i,del3_W_min,del3_W_max,del3_W_L, del3_W_R,del3_W_C,d,a_dx);
			// }
			
			
			// Slope flattening starts here
			Scalar eta ;
			Scalar eta_old ;

			for (int d = 0; d < DIM; d++)
			{

				//Vector W_ave_ahead = alias(a_W_ave,Point::Basis(d)*(-1));
				Vector W_ave_ahead = m_ahead_shift[d](a_W_ave);
				//Vector W_ave_ahead2 = alias(a_W_ave,Point::Basis(d)*(-2));
				Vector W_ave_ahead2 = m_ahead_shift[d](W_ave_ahead);
				//Vector W_ave_ahead3 = alias(a_W_ave,Point::Basis(d)*(-3));
				Vector W_ave_ahead3 = m_ahead_shift[d](W_ave_ahead2);
				
				//Vector W_ave_behind = alias(a_W_ave, Point::Basis(d)*(1));
				Vector W_ave_behind = m_behind_shift[d](a_W_ave);
				//Vector W_ave_behind2 = alias(a_W_ave,Point::Basis(d)*(2));	
				Vector W_ave_behind2 = m_behind_shift[d](W_ave_behind);	
				//Vector W_ave_behind3 = alias(a_W_ave,Point::Basis(d)*(3));
				Vector W_ave_behind3 = m_behind_shift[d](W_ave_behind2);	
				
				Scalar eta_tilde_d = forall<double>( eta_tilde_d_calc, a_W_ave, W_ave_ahead, W_ave_ahead2,
											  W_ave_behind, W_ave_behind2, d);
										  
				//Scalar eta_tilde_d_ahead = alias(eta_tilde_d,Point::Basis(d)*(-1));
				Scalar eta_tilde_d_ahead = m_ahead_shift[d](eta_tilde_d);
				//Scalar eta_tilde_d_behind = alias(eta_tilde_d,Point::Basis(d)*(1));	
				Scalar eta_tilde_d_behind = m_behind_shift[d](eta_tilde_d);	
				
				Scalar eta_d = forall<double>( eta_d_calc, eta_tilde_d, eta_tilde_d_ahead, eta_tilde_d_behind);
				if (d>0){
					eta = forall<double>(eta_calc, eta_d, eta_old);
				} else {
					eta = m_copy(eta_d);
				}	
				eta_old = m_copy(eta);
			}
			
			
			Vector W_ave_low_ahead_lim_flat = forall<double,NUMCOMPS>(Flat_calc, W_ave_low_ahead_limited, a_W_ave, eta, a_d);
			a_W_ave_high_lim_flat = forall<double,NUMCOMPS>(Flat_calc, W_ave_high_limited, a_W_ave, eta, a_d);
			a_W_ave_low_lim_flat = m_behind_shift[a_d](W_ave_low_ahead_limited);
			//Slope flattening ends here
	}
	
}