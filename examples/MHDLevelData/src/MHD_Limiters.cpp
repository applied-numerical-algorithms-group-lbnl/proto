#include "Proto.H"
#include "MHD_Limiters.H"
#include "CommonTemplates.H"
#include "Proto_Timer.H"
#include "Proto_WriteBoxData.H"
#include "MHD_Input_Parsing.H"
#include <iomanip>

typedef BoxData<double,1,HOST> Scalar;
typedef BoxData<double,NUMCOMPS,HOST> Vector;

extern Parsefrominputs inputs;

namespace MHD_Limiters {

	PROTO_KERNEL_START
	void del2_W_c_calcF(State& del2_W_c,
	                    const State& a_W_ave,
	                    const State& a_W_ave_behind,
	                    const State& a_W_ave_ahead)
	{
		for (int i=0; i< NUMCOMPS; i++) {
			del2_W_c(i) = a_W_ave_behind(i) - 2.0*a_W_ave(i) + a_W_ave_ahead(i);
		}
	}
	PROTO_KERNEL_END(del2_W_c_calcF, del2_W_c_calc)


	PROTO_KERNEL_START
	void del3_W_calcF(State& a_del3_W,
	                  const State& a_del2_W_c,
	                  const State& a_del2_W_c_behind)
	{
		for (int i=0; i< NUMCOMPS; i++) {
			a_del3_W(i) = a_del2_W_c(i) - a_del2_W_c_behind(i);
		}
	}
	PROTO_KERNEL_END(del3_W_calcF, del3_W_calc)



	const bool printlim = false;

	PROTO_KERNEL_START
	void limiter_calcF(
		const Point& a_pt,
		State& a_W_ave_low_ahead_limited,
		State& a_W_ave_high_limited,
		const State& a_W_ave_low_ahead,
		const State& a_W_ave_high,
		const State& a_W_ave,
		const State& a_W_ave_ahead,
		const State& a_W_ave_ahead2,
		const State& a_W_ave_behind,
		const State& a_W_ave_behind2,
		const State& a_del2_W_c,
		const State& a_del2_W_c_ahead,
		const State& a_del2_W_c_behind,
		const State& a_del3_W,
		const State& a_del3_W_L,
		const State& a_del3_W_R,
		const State& a_del3_W_behind,
		const State& a_del3_W_ahead,
		const State& a_del3_W_ahead2,
		int a_dir,
		const double a_dx,
		const double a_dy,
		const double a_dz)
	{
		//if (a_pt [1] < 0) cout << a_pt[0] << " "<< a_pt [1] << endl;
		double rhs = 1.0 - 1.0e-12;
		double rhs_test, lhs_test, rhs_test2, lhs_test2,rhs_test_a, lhs_test_a, rhs_rho, lhs_rho;
		double a_del_W_f_m, a_del_W_f_p, a_del2_W_f, a_del3_W_C, a_del2_W_lim, a_rho_i, a_del3_W_min, a_del3_W_max;
		for (int i=0; i< NUMCOMPS; i++)
		{
		//if (i == 3){
			a_del_W_f_m = a_W_ave(i) - a_W_ave_high(i);
			a_del_W_f_p = a_W_ave_low_ahead(i) - a_W_ave(i);
			a_del2_W_f = 6.0*(a_W_ave_high(i) - 2.0*a_W_ave(i) + a_W_ave_low_ahead(i));
			a_del3_W_C = (a_del2_W_c_ahead(i) - a_del2_W_c_behind(i))/2.0;
			if ((a_del2_W_c_behind(i) >= 0.0) && (a_del2_W_c(i) >= 0.0) && (a_del2_W_c_ahead(i) >= 0.0) && (a_del2_W_f >= 0.0)) {
				a_del2_W_lim = std::min({std::abs(a_del2_W_f), 1.25*std::abs(a_del2_W_c_behind(i)),
				                         1.25*std::abs(a_del2_W_c(i)), 1.25*std::abs(a_del2_W_c_ahead(i))});
			} else if ((a_del2_W_c_behind(i) < 0.0) && (a_del2_W_c(i) < 0.0) && (a_del2_W_c_ahead(i) < 0.0) && (a_del2_W_f < 0.0)) {
				a_del2_W_lim = -1.0*std::min({std::abs(a_del2_W_f), 1.25*std::abs(a_del2_W_c_behind(i)),
				                              1.25*std::abs(a_del2_W_c(i)), 1.25*std::abs(a_del2_W_c_ahead(i))});
			} else {
				a_del2_W_lim = 0.0;
			}
			rhs_rho = 1.0e-12*std::max({std::abs(a_W_ave_behind2(i)), std::abs(a_W_ave_behind(i)), std::abs(a_W_ave(i)), std::abs(a_W_ave_ahead(i)), std::abs(a_W_ave_ahead2(i))});
			lhs_rho = std::abs(a_del2_W_f);
			if (lhs_rho <= rhs_rho) {
				a_rho_i = 0.0;
			} else {
				a_rho_i = a_del2_W_lim/a_del2_W_f;
			}
			a_del3_W_min = std::min({a_del3_W_behind(i), a_del3_W(i), a_del3_W_ahead(i), a_del3_W_ahead2(i)});
			a_del3_W_max = std::max({a_del3_W_behind(i), a_del3_W(i), a_del3_W_ahead(i), a_del3_W_ahead2(i)});



			if (printlim)
			{
				if (i==1 && a_dir == 0){
					if (a_pt[0]  >= 16 && a_pt[0]  <= 25 && a_pt [1] == 45 ) {
						cout << setw(30) << setprecision(18) << a_pt[0]
						     << setw(30) << setprecision(18) << a_pt[1]
							 << setw(30) << setprecision(18) << a_W_ave(i)
						     << setw(30) << setprecision(18) << a_W_ave_low_ahead(i)
						     << setw(30) << setprecision(18) << a_W_ave_high(i)
						     << setw(30) << setprecision(18) << a_del_W_f_p
						     << setw(30) << setprecision(18) << a_del_W_f_m
						     << setw(30) << setprecision(18) << a_del2_W_c(i)
						     << setw(30) << setprecision(18) << a_del2_W_f
						     << setw(30) << setprecision(18) << a_del3_W(i)
						     << setw(30) << setprecision(18) << a_del2_W_lim
						     << setw(30) << setprecision(18) << a_rho_i
						     << setw(30) << setprecision(18) << a_del3_W_min
						     << setw(30) << setprecision(18) << a_del3_W_max
						     << setw(30) << setprecision(18) << a_del3_W_L(i)
						     << setw(30) << setprecision(18) << a_del3_W_R(i)
						     << setw(30) << setprecision(18) << a_del3_W_C
							 << endl;
					}
				}
			}




			if ((a_del_W_f_m * a_del_W_f_p <= 0.0) || ((a_W_ave(i) - a_W_ave_behind2(i))*(a_W_ave_ahead2(i) - a_W_ave(i)) <= 0.0)) {
				if (a_rho_i < rhs) {
					rhs_test = a_del3_W_max - a_del3_W_min;
					lhs_test = 0.1*std::max({std::abs(a_del3_W_min), std::abs(a_del3_W_max)});
					lhs_test_a = std::min({std::abs(a_del3_W_L(i)),std::abs(a_del3_W_R(i))})*1.1;
					rhs_test_a = std::abs(a_del3_W_C);
					if (lhs_test <= rhs_test){
					//if (!(lhs_test_a > rhs_test_a && a_del3_W_L(i)*a_del3_W_R(i) > 0.)) {

						rhs_test2 = 2.0*std::abs(a_del_W_f_m);
						lhs_test2 = std::abs(a_del_W_f_p);
						if (a_del_W_f_m * a_del_W_f_p < 0.0) {
							a_W_ave_low_ahead_limited(i) = a_W_ave(i) + a_rho_i*a_del_W_f_p;

							// if (std::abs(a_W_ave_low_ahead(i) - a_W_ave_low_ahead_limited(i)) > 1.0e-15){
							// // if (a_W_ave_low_ahead(i) != a_W_ave_low_ahead_limited(i)){
							// 	// cout << "Lo Reason 1 " << "i = " << a_pt[0] << " j = " << a_pt[1] << " Comp = " << i << " Diff = " << setw(30) << setprecision(18) << std::abs(a_W_ave_low_ahead(i) - a_W_ave_low_ahead_limited(i)) << setw(30) << setprecision(18) << a_W_ave_low_ahead(i) << setw(30) << setprecision(18) << a_W_ave_low_ahead_limited(i) << endl;
							// 	// a_W_ave_low_ahead_limited(i) = a_W_ave_low_ahead(i);
							// }
						} else if (lhs_test2 >= rhs_test2) {
							a_W_ave_low_ahead_limited(i) = a_W_ave(i) + 2.0*(1.0-a_rho_i)*a_del_W_f_m + a_rho_i*a_del_W_f_p;

							// if (std::abs(a_W_ave_low_ahead(i) - a_W_ave_low_ahead_limited(i)) > 1.0e-15){
							// // if (a_W_ave_low_ahead(i) != a_W_ave_low_ahead_limited(i)){
							// 	// cout << "Lo Reason 1 " << "i = " << a_pt[0] << " j = " << a_pt[1] << " Comp = " << i << " Diff = " << setw(30) << setprecision(18) << std::abs(a_W_ave_low_ahead(i) - a_W_ave_low_ahead_limited(i)) << setw(30) << setprecision(18) << a_W_ave_low_ahead(i) << setw(30) << setprecision(18) << a_W_ave_low_ahead_limited(i) << endl;
							// 	// a_W_ave_low_ahead_limited(i) = a_W_ave_low_ahead(i);
							// }
						}

						rhs_test2 = 2.0*std::abs(a_del_W_f_p);
						lhs_test2 = std::abs(a_del_W_f_m);
						if (a_del_W_f_m * a_del_W_f_p < 0.0) {
							a_W_ave_high_limited(i) = a_W_ave(i) - a_rho_i*a_del_W_f_m;

							// if (std::abs(a_W_ave_high(i) - a_W_ave_high_limited(i)) > 1.0e-15){
							// 	// cout << "Hi Reason 1 " << "i = " << a_pt[0] << " j = " << a_pt[1] << " Comp = " << i << " Diff = " << setw(30) << setprecision(18) << std::abs(a_W_ave_high(i) - a_W_ave_high_limited(i)) << setw(30) << setprecision(18) << a_W_ave_high(i) << setw(30) << setprecision(18) << a_W_ave_high_limited(i) << endl;
							// 	// a_W_ave_high_limited(i) = a_W_ave_high(i);
							// }
						} else if (lhs_test2 >= rhs_test2) {
							a_W_ave_high_limited(i) = a_W_ave(i) - 2.0*(1.0-a_rho_i)*a_del_W_f_p - a_rho_i*a_del_W_f_m;

							// if (std::abs(a_W_ave_high(i) - a_W_ave_high_limited(i)) > 1.0e-15){
							// 	// cout << "Hi Reason 2 " << "i = " << a_pt[0] << " j = " << a_pt[1] << " Comp = " << i << " Diff = " << setw(30) << setprecision(18) << std::abs(a_W_ave_high(i) - a_W_ave_high_limited(i)) << setw(30) << setprecision(18) << a_W_ave_high(i) << setw(30) << setprecision(18) << a_W_ave_high_limited(i) << endl;
							// 	// a_W_ave_high_limited(i) = a_W_ave_high(i);
							// }
						}

					}
				}
			} else {

				rhs_test2 = 2.0*std::abs(a_del_W_f_m);
				lhs_test2 = std::abs(a_del_W_f_p);
				if (lhs_test2 >= rhs_test2) {
					a_W_ave_low_ahead_limited(i) = a_W_ave(i) + 2.0*a_del_W_f_m;



					// if (std::abs(a_W_ave_low_ahead(i) - a_W_ave_low_ahead_limited(i)) > 1.0e-15){
					// // if (a_W_ave_low_ahead(i) != a_W_ave_low_ahead_limited(i)){
					// 	// cout << "Lo Reason 1 " << "i = " << a_pt[0] << " j = " << a_pt[1] << " Comp = " << i << " Diff = " << setw(30) << setprecision(18) << std::abs(a_W_ave_low_ahead(i) - a_W_ave_low_ahead_limited(i)) << setw(30) << setprecision(18) << a_W_ave_low_ahead(i) << setw(30) << setprecision(18) << a_W_ave_low_ahead_limited(i) << endl;
					// 	// a_W_ave_low_ahead_limited(i) = a_W_ave_low_ahead(i);
					// }
				}

				rhs_test2 = 2.0*std::abs(a_del_W_f_p);
				lhs_test2 = std::abs(a_del_W_f_m);
				if (lhs_test2 >= rhs_test2) {
					a_W_ave_high_limited(i) = a_W_ave(i) - 2.0*a_del_W_f_p;



					// if (std::abs(a_W_ave_high(i) - a_W_ave_high_limited(i)) > 1.0e-15){
					// 	// cout << "Hi Reason 3 " << "i = " << a_pt[0] << " j = " << a_pt[1] << " Comp = " << i  << " Diff = " << setw(30) << setprecision(18) << std::abs(a_W_ave_high(i) - a_W_ave_high_limited(i)) << setw(30) << setprecision(18) << a_W_ave_high(i) << setw(30) << setprecision(18) << a_W_ave_high_limited(i) << endl;
					// 	// a_W_ave_high_limited(i) = a_W_ave_high(i);
					// }
				}

			}
		//}

		}

	}
	PROTO_KERNEL_END(limiter_calcF, limiter_calc)



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
		if (arg > z1) {si = 0.0;}
		if (arg < z0) {si = 1.0;}
		if ((arg<=z1) && (arg>=z0)) {si = 1.0 - (arg-z0)/(z1-z0);}

		double lhs1 = a_W_ave_behind(1+a_d)-a_W_ave_ahead(1+a_d);
		double lhs2 = std::abs(p_ahead-p_behind)/std::min({p_ahead,p_behind});
		if ((lhs1 > 0.0) && (lhs2 > delta)) {
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
	                int a_dir)
	{

		for (int i=0; i< NUMCOMPS; i++) {
			//if (i != (NUMCOMPS + a_dir - 2))
			if (i != 100)
			{
				a_flattened(i) = (1.0-a_eta(0))*a_not_flattened(i) + a_eta(0)*a_W_ave(i);       //This works
			} else {
				a_flattened(i) = a_not_flattened(i);
			}
		}

	}
	PROTO_KERNEL_END(Flat_calcF, Flat_calc)


	PROTO_KERNEL_START
	void compare_calcF(const Point& a_pt,
					State& a_W_diff,
	                const State& a_W_ave_high_limited,
	                const State& a_W_ave)
	{

		for (int i=0; i< NUMCOMPS; i++) {
			if (a_W_ave_high_limited(i) != a_W_ave(i)) cout << a_pt[0] << a_pt[1] << endl;
		}

	}
	PROTO_KERNEL_END(compare_calcF, compare_calc)

	PROTO_KERNEL_START
	void ResetpolesF(const Point& a_pt,
					State& a_W_lim,
	                const State& a_W)
	{

		for (int i=0; i< NUMCOMPS; i++) {
			a_W_lim(i) = a_W(i);
		}

	}
	PROTO_KERNEL_END(ResetpolesF, Resetpoles)


	void MHD_Limiters(BoxData<double,NUMCOMPS>& a_W_ave_low_lim_flat,
	                  BoxData<double,NUMCOMPS>& a_W_ave_high_lim_flat,
	                  BoxData<double,NUMCOMPS>& a_W_ave_low,
	                  BoxData<double,NUMCOMPS>& a_W_ave_high,
	                  BoxData<double,NUMCOMPS>& a_W_ave,
	                  BoxData<double,NUMCOMPS>& a_W_bar,
	                  const int a_d,
	                  const double a_dx,
	                  const double a_dy,
	                  const double a_dz)
	{

		static Stencil<double> m_copy;
		m_copy = 1.0*Shift(Point::Zeros());

		//Limiter Starts here
		Vector W_ave_low_ahead = alias(a_W_ave_low,Point::Basis(a_d)*(-1));
		//Vector W_ave_low_ahead_limited = alias(W_ave_low_ahead);  // This works, but won't be able to compare diffrences as if W_ave_low_ahead_limited changes, so does W_ave_low_ahead
		Vector W_ave_low_ahead_limited = m_copy(W_ave_low_ahead);
		//Vector W_ave_high_limited = alias(a_W_ave_high);
		Vector W_ave_high_limited = m_copy(a_W_ave_high);
		if (inputs.limiter_apply == 1) {

			Vector W_ave_ahead = alias(a_W_ave,Point::Basis(a_d)*(-1));
			Vector W_ave_ahead2 = alias(a_W_ave,Point::Basis(a_d)*(-2));
			Vector W_ave_behind = alias(a_W_ave,Point::Basis(a_d)*(1));
			Vector W_ave_behind2 = alias(a_W_ave,Point::Basis(a_d)*(2));
			Vector del2_W_c  = forall<double,NUMCOMPS>( del2_W_c_calc, a_W_ave, W_ave_behind, W_ave_ahead);
			Vector del2_W_c_ahead = alias(del2_W_c,Point::Basis(a_d)*(-1));
			Vector del2_W_c_behind = alias(del2_W_c,Point::Basis(a_d)*(1));
			Vector del3_W = forall<double,NUMCOMPS>( del3_W_calc, del2_W_c, del2_W_c_behind);
			Vector del3_W_L = forall<double,NUMCOMPS>( del3_W_calc, del2_W_c, del2_W_c_behind);
			Vector del3_W_R = forall<double,NUMCOMPS>( del3_W_calc, del2_W_c_ahead, del2_W_c);
			Vector del3_W_behind = alias(del3_W,Point::Basis(a_d)*(1));
			Vector del3_W_ahead = alias(del3_W,Point::Basis(a_d)*(-1));
			Vector del3_W_ahead2 = alias(del3_W,Point::Basis(a_d)*(-2));


		forallInPlace_p( limiter_calc,W_ave_low_ahead_limited, W_ave_high_limited, W_ave_low_ahead, a_W_ave_high,
			                 a_W_ave,W_ave_ahead, W_ave_ahead2, W_ave_behind, W_ave_behind2, del2_W_c, del2_W_c_ahead, del2_W_c_behind, del3_W,
			                 del3_W_L, del3_W_R, del3_W_behind, del3_W_ahead, del3_W_ahead2, a_d, a_dx, a_dy, a_dz);
		}
		Box dbx0 = W_ave_low_ahead_limited.box();
		Box dbx1 = W_ave_low_ahead.box();

		//cout << "Original box: " << dbx1.low()[0] << " " << dbx1.high()[0] << " "  << dbx1.low()[1] << " "  << dbx1.high()[1] << endl;
		//cout << "Limited box: " << dbx0.low()[0] << " "  << dbx0.high()[0] << " "  << dbx0.low()[1] << " "  << dbx0.high()[1] << endl;
		Vector W_diff(dbx0);
		//forallInPlace_p(compare_calc, W_diff, W_ave_high_limited, a_W_ave_high);

		// Slope flattening starts here
		if (inputs.slope_flattening_apply == 1) {
			Scalar eta;
			Scalar eta_old;

			for (int d = 0; d < DIM; d++)
			{

				Vector W_bar_ahead = alias(a_W_bar,Point::Basis(d)*(-1));
				Vector W_bar_ahead2 = alias(a_W_bar,Point::Basis(d)*(-2));
				Vector W_bar_behind = alias(a_W_bar,Point::Basis(d)*(1));
				Vector W_bar_behind2 = alias(a_W_bar,Point::Basis(d)*(2));

				Scalar eta_tilde_d = forall<double>( eta_tilde_d_calc, a_W_bar, W_bar_ahead, W_bar_ahead2,
				                                     W_bar_behind, W_bar_behind2, d);

				Scalar eta_tilde_d_ahead = alias(eta_tilde_d,Point::Basis(d)*(-1));
				Scalar eta_tilde_d_behind = alias(eta_tilde_d,Point::Basis(d)*(1));

				Scalar eta_d = forall<double>( eta_d_calc, eta_tilde_d, eta_tilde_d_ahead, eta_tilde_d_behind);
				if (d>0) {
					eta = forall<double>(eta_calc, eta_d, eta_old);
				} else {
					eta = alias(eta_d);
				}
				eta_old = alias(eta);
			}
			Vector W_ave_low_ahead_lim_flat = forall<double,NUMCOMPS>(Flat_calc, W_ave_low_ahead_limited, a_W_ave, eta, a_d);
			a_W_ave_high_lim_flat = forall<double,NUMCOMPS>(Flat_calc, W_ave_high_limited, a_W_ave, eta, a_d);
			a_W_ave_low_lim_flat = alias(W_ave_low_ahead_lim_flat,Point::Basis(a_d)*(1));
			//Slope flattening ends here
		} else {
			Vector W_ave_low_ahead_lim_flat = alias(W_ave_low_ahead_limited);
			a_W_ave_high_lim_flat = alias(W_ave_high_limited);
			a_W_ave_low_lim_flat = alias(W_ave_low_ahead_lim_flat,Point::Basis(a_d)*(1));

#if DIM == 3
			// if (a_d == 1){
			// if (inputs.grid_type_global == 2){
			// 	Box bb = a_W_ave_high_lim_flat.box();
			// 	if (bb.low()[1] <= 0){
			// 		Point pole_lo = Point(bb.low()[0], 0, bb.low()[2]);
			// 		Point pole_hi = Point(bb.high()[0], 0, bb.high()[2]);
			// 		Box BoundBox(pole_lo,pole_hi);
			// 		forallInPlace_p(Resetpoles,BoundBox,a_W_ave_low_lim_flat,a_W_ave_low);
			// 	}

			// 	bb = a_W_ave_high_lim_flat.box();
			// 	if (bb.high()[1] >= 63){
			// 		Point pole_lo = Point(bb.low()[0], 63, bb.low()[2]);
			// 		Point pole_hi = Point(bb.high()[0], 63, bb.high()[2]);
			// 		Box BoundBox(pole_lo,pole_hi);
			// 		forallInPlace_p(Resetpoles,BoundBox,a_W_ave_high_lim_flat,a_W_ave_high);
			// 	}
			// }
			// }
#endif
		}


	}

}
