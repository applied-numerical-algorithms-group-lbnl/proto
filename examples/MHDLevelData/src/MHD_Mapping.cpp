#include "Proto.H"
#include "MHD_Mapping.H"
#include "CommonTemplates.H"
#include "Proto_Timer.H"
#include "Proto_WriteBoxData.H"
#include "MHDOp.H"

const double C1_fix = 0.1; // A constant in wavy grid definition.
const double r_in = 0.2;
const double r_out = 0.8;
const double C_rad = 1.0; // A constant in exponential dr in spherical grid.
extern int grid_type_global;
typedef BoxData<double,1,1,1> Scalar;
typedef BoxData<double,NUMCOMPS,1,1> Vector;

namespace MHD_Mapping {

	PROTO_KERNEL_START
	void iotaFuncF(Point           & a_p,
	               V               & a_X,
	               const double a_dx,
	               const double a_dy,
	               const double a_dz)
	{
		for (int ii = 0; ii < DIM; ii++)
		{
			double dxd;
			if (ii == 0) dxd = a_dx;
			if (ii == 1) dxd = a_dy;
			if (ii == 2) dxd = a_dz;
			a_X(ii) = a_p[ii]*dxd + 0.5*dxd;
		}
	}
	PROTO_KERNEL_END(iotaFuncF,iotaFunc)


	void eta_calc(BoxData<double,DIM>& a_eta,
	              const Box& a_bx,
	              const double a_dx,
	              const double a_dy,
	              const double a_dz)
	{
		forallInPlace_p(iotaFunc, a_bx, a_eta, a_dx, a_dy, a_dz);
	}


	PROTO_KERNEL_START
	void eta_to_xF(Var<double,DIM>& a_x,
	               const Var<double,DIM>& a_eta)
	{

#if DIM == 2
		if (grid_type_global == 0) {
			a_x(0) = a_eta(0);
			a_x(1) = a_eta(1);
		}
		if (grid_type_global == 1) {
			double C1 = C1_fix;
			double C2 = C1;
			a_x(0) = a_eta(0) + C1*sin(2.0*PI*a_eta(0))*sin(2.0*PI*a_eta(1));
			a_x(1) = a_eta(1) + C2*sin(2.0*PI*a_eta(0))*sin(2.0*PI*a_eta(1));
		}
		if (grid_type_global == 2) {
			a_x(0) = (r_in + a_eta(0)*(r_out-r_in))*cos(2.0*PI*a_eta(1));
			a_x(1) = (r_in + a_eta(0)*(r_out-r_in))*sin(2.0*PI*a_eta(1));
		}
#endif
#if DIM == 3
		if (grid_type_global == 0) {
			a_x(0) = a_eta(0);
			a_x(1) = a_eta(1);
			a_x(2) = a_eta(2);
		}
		if (grid_type_global == 2){
			double R_t = (r_out - r_in)/(exp(C_rad) - 1.0);
			double r = r_in + R_t*(exp(C_rad*a_eta(0)) - 1.0);
			a_x(0) = r*sin(PI*a_eta(2))*cos(2.0*PI*a_eta(1));
			a_x(1) = r*sin(PI*a_eta(2))*sin(2.0*PI*a_eta(1));
			a_x(2) = r*cos(PI*a_eta(2));
		}

#endif
	}
	PROTO_KERNEL_END(eta_to_xF, eta_to_x)


	void eta_to_x_calc(BoxData<double,DIM>& a_x,
	                   const BoxData<double,DIM>& a_eta)
	{
		a_x = forall<double,DIM>(eta_to_x,a_eta);
	}

	void phys_coords_calc(BoxData<double,DIM>& a_x,
	                      const Box& dbx1,
	                      const double a_dx,
	                      const double a_dy,
	                      const double a_dz)
	{
		Box dbx0=a_x.box();
		BoxData<double,DIM> eta(dbx0);
		MHD_Mapping::eta_calc(eta,dbx0,a_dx, a_dy, a_dz);
		MHD_Mapping::eta_to_x_calc(a_x,eta);
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
	                    const double a_dx,
	                    const double a_dy,
	                    const double a_dz)
	{

#if DIM == 2
		double x_hi = (a_pt[0] + 1.0)*a_dx;
		double x_lo = (a_pt[0] - 0.0)*a_dx;
		double y_hi = (a_pt[1] + 1.0)*a_dy;
		double y_lo = (a_pt[1] - 0.0)*a_dy;

		if (grid_type_global == 0) {
			if (a_s == 0 && a_d == 0) a_X_ave_f(0) = ((x_lo*y_hi)
				                                  -(x_lo*y_lo))/DIM/a_dy;
			if (a_s == 1 && a_d == 0) a_X_ave_f(0) = ((y_hi*y_hi/2.0)
				                                  -(y_lo*y_lo/2.0))/DIM/a_dy;
			if (a_s == 0 && a_d == 1) a_X_ave_f(0) = ((x_hi*x_hi/2.0)
				                                  -(x_lo*x_lo/2.0))/DIM/a_dx;
			if (a_s == 1 && a_d == 1) a_X_ave_f(0) = ((x_hi*y_lo)
				                                  -(x_lo*y_lo))/DIM/a_dx;
		}

		if (grid_type_global == 1) {
			double C1 = C1_fix;
			double C2 = C1;
			if (a_s == 0 && a_d == 0) a_X_ave_f(0) = ((x_lo*y_hi-(C1/2./PI)*sin(2.*PI*x_lo)*cos(2.*PI*y_hi))
				                                  -(x_lo*y_lo-(C1/2./PI)*sin(2.*PI*x_lo)*cos(2.*PI*y_lo)))/DIM/a_dy;
			if (a_s == 1 && a_d == 0) a_X_ave_f(0) = ((y_hi*y_hi/2.0-(C2/2./PI)*sin(2.*PI*x_lo)*cos(2.*PI*y_hi))
				                                  -(y_lo*y_lo/2.0-(C2/2./PI)*sin(2.*PI*x_lo)*cos(2.*PI*y_lo)))/DIM/a_dy;
			if (a_s == 0 && a_d == 1) a_X_ave_f(0) = ((x_hi*x_hi/2.0-(C1/2./PI)*cos(2.*PI*x_hi)*sin(2.*PI*y_lo))
				                                  -(x_lo*x_lo/2.0-(C1/2./PI)*cos(2.*PI*x_lo)*sin(2.*PI*y_lo)))/DIM/a_dx;
			if (a_s == 1 && a_d == 1) a_X_ave_f(0) = ((x_hi*y_lo-(C2/2./PI)*cos(2.*PI*x_hi)*sin(2.*PI*y_lo))
				                                  -(x_lo*y_lo-(C2/2./PI)*cos(2.*PI*x_lo)*sin(2.*PI*y_lo)))/DIM/a_dx;
		}

		if (grid_type_global == 2) {
			double r_diff = r_out-r_in;
			if (a_s == 0 && a_d == 0) a_X_ave_f(0) = (((r_in+x_lo*r_diff)*sin(2.*PI*y_hi)/2./PI)
				                                  -((r_in+x_lo*r_diff)*sin(2.*PI*y_lo)/2./PI))/DIM/a_dy;
			if (a_s == 1 && a_d == 0) a_X_ave_f(0) = -(((r_in+x_lo*r_diff)*cos(2.*PI*y_hi)/2./PI)
				                                   -((r_in+x_lo*r_diff)*cos(2.*PI*y_lo)/2./PI))/DIM/a_dy;
			if (a_s == 0 && a_d == 1) a_X_ave_f(0) = (((r_in*x_hi+x_hi*x_hi*r_diff/2.0)*cos(2.*PI*y_lo))
				                                  -((r_in*x_lo+x_lo*x_lo*r_diff/2.0)*cos(2.*PI*y_lo)))/DIM/a_dx;
			if (a_s == 1 && a_d == 1) a_X_ave_f(0) = (((r_in*x_hi+x_hi*x_hi*r_diff/2.0)*sin(2.*PI*y_lo))
				                                  -((r_in*x_lo+x_lo*x_lo*r_diff/2.0)*sin(2.*PI*y_lo)))/DIM/a_dx;
		}
#endif

#if DIM == 3


		if (grid_type_global == 0) {
			if (a_s == 0 && a_d == 0) a_X_ave_f(0) = (a_pt[0] + 0.0)*a_dx/DIM;
			if (a_s == 1 && a_d == 0) a_X_ave_f(0) = (a_pt[1] + 0.5)*a_dy/DIM;
			if (a_s == 2 && a_d == 0) a_X_ave_f(0) = (a_pt[2] + 0.5)*a_dz/DIM;
			if (a_s == 0 && a_d == 1) a_X_ave_f(0) = (a_pt[0] + 0.5)*a_dx/DIM;
			if (a_s == 1 && a_d == 1) a_X_ave_f(0) = (a_pt[1] + 0.0)*a_dy/DIM;
			if (a_s == 2 && a_d == 1) a_X_ave_f(0) = (a_pt[2] + 0.5)*a_dz/DIM;
			if (a_s == 0 && a_d == 2) a_X_ave_f(0) = (a_pt[0] + 0.5)*a_dx/DIM;
			if (a_s == 1 && a_d == 2) a_X_ave_f(0) = (a_pt[1] + 0.5)*a_dy/DIM;
			if (a_s == 2 && a_d == 2) a_X_ave_f(0) = (a_pt[2] + 0.0)*a_dz/DIM;
		}

		if (grid_type_global == 2) {
			double eta1_hi = (a_pt[0] + 1.0)*a_dx;
			double eta1_lo = (a_pt[0] - 0.0)*a_dx;
			double eta2_hi = (a_pt[1] + 1.0)*a_dy;
			double eta2_lo = (a_pt[1] - 0.0)*a_dy;
			double eta3_hi = (a_pt[2] + 1.0)*a_dz;
			double eta3_lo = (a_pt[2] - 0.0)*a_dz;

			double R_t = (r_out - r_in)/(exp(C_rad) - 1.0);
			double r_lo = r_in + R_t*(exp(C_rad*eta1_lo) - 1.0);
			double coseta3_hi = cos(PI*eta3_hi);
			double coseta3_lo = cos(PI*eta3_lo);
			double coseta2_hi = cos(2.0*PI*eta2_hi);
			double coseta2_lo = cos(2.0*PI*eta2_lo);
			double sineta3_hi = sin(PI*eta3_hi);
			double sineta3_lo = sin(PI*eta3_lo);
			double sineta2_hi = sin(2.0*PI*eta2_hi);
			double sineta2_lo = sin(2.0*PI*eta2_lo);

			double r_int_lo = r_in*eta1_lo + R_t*((exp(C_rad*eta1_lo)/C_rad) - eta1_lo);
			double r_int_hi = r_in*eta1_hi + R_t*((exp(C_rad*eta1_hi)/C_rad) - eta1_hi);

			if (a_s == 0 && a_d == 0) a_X_ave_f(0) = r_lo*(coseta3_lo-coseta3_hi)*(sineta2_hi-sineta2_lo)/DIM/2.0/PI/PI/a_dy/a_dz;
			if (a_s == 1 && a_d == 0) a_X_ave_f(0) = r_lo*(coseta3_lo-coseta3_hi)*(coseta2_lo-coseta2_hi)/DIM/2.0/PI/PI/a_dy/a_dz;
			if (a_s == 2 && a_d == 0) a_X_ave_f(0) = r_lo*(sineta3_hi-sineta3_lo)/DIM/PI/a_dz;
			if (a_s == 0 && a_d == 1) a_X_ave_f(0) = coseta2_lo*(r_int_hi-r_int_lo)*(coseta3_lo-coseta3_hi)/DIM/PI/a_dx/a_dz;
			if (a_s == 1 && a_d == 1) a_X_ave_f(0) = sineta2_lo*(r_int_hi-r_int_lo)*(coseta3_lo-coseta3_hi)/DIM/PI/a_dx/a_dz;
			if (a_s == 2 && a_d == 1) a_X_ave_f(0) = (r_int_hi-r_int_lo)*(sineta3_hi-sineta3_lo)/DIM/PI/a_dx/a_dz;
			if (a_s == 0 && a_d == 2) a_X_ave_f(0) = sineta3_lo*(r_int_hi-r_int_lo)*(sineta2_hi-sineta2_lo)/DIM/2.0/PI/a_dx/a_dy;
			if (a_s == 1 && a_d == 2) a_X_ave_f(0) = sineta3_lo*(r_int_hi-r_int_lo)*(coseta2_lo-coseta2_hi)/DIM/2.0/PI/a_dx/a_dy;
			if (a_s == 2 && a_d == 2) a_X_ave_f(0) = coseta3_lo*(r_int_hi-r_int_lo)/DIM/a_dx;
		}

#endif

	}
	PROTO_KERNEL_END(X_ave_f_calcF, X_ave_f_calc)



	PROTO_KERNEL_START
	void N_ave_f_calcF( const Point& a_pt,
	                    Var<double,1>& a_N_ave_f,
	                    int a_s,
	                    int a_d,
	                    const double a_dx,
	                    const double a_dy,
	                    const double a_dz)
	{
#if DIM == 2
		int x_loc = a_pt[0];
		int y_loc = a_pt[1];
		double x_hi = (x_loc + 1.0)*a_dx;
		double x_lo = (x_loc - 0.0)*a_dx;
		double y_hi = (y_loc + 1.0)*a_dy;
		double y_lo = (y_loc - 0.0)*a_dy;

		if (grid_type_global == 0) {
			if (a_s == 0 && a_d == 0) a_N_ave_f(0) = 1.0;
			if (a_s == 1 && a_d == 0) a_N_ave_f(0) = 0.0;
			if (a_s == 0 && a_d == 1) a_N_ave_f(0) = 0.0;
			if (a_s == 1 && a_d == 1) a_N_ave_f(0) = 1.0;
		}

		if (grid_type_global == 1) {
			double C1 = C1_fix;
			double C2 = C1;
			if (a_s == 0 && a_d == 0) a_N_ave_f(0) = ((y_hi+C2*sin(2.*PI*x_lo)*sin(2.*PI*y_hi))
				                                  -(y_lo+C2*sin(2.*PI*x_lo)*sin(2.*PI*y_lo)))/a_dy;
			if (a_s == 1 && a_d == 0) a_N_ave_f(0) = -((x_lo+C1*sin(2.*PI*x_lo)*sin(2.*PI*y_hi))
				                                   -(x_lo+C1*sin(2.*PI*x_lo)*sin(2.*PI*y_lo)))/a_dy;
			if (a_s == 0 && a_d == 1) a_N_ave_f(0) = -((y_lo+C2*sin(2.*PI*x_hi)*sin(2.*PI*y_lo))
				                                   -(y_lo+C2*sin(2.*PI*x_lo)*sin(2.*PI*y_lo)))/a_dx;
			if (a_s == 1 && a_d == 1) a_N_ave_f(0) = ((x_hi+C1*sin(2.*PI*x_hi)*sin(2.*PI*y_lo))
				                                  -(x_lo+C1*sin(2.*PI*x_lo)*sin(2.*PI*y_lo)))/a_dx;
		}

		if (grid_type_global == 2) {
			double r_diff = r_out-r_in;
			if (a_s == 0 && a_d == 0) a_N_ave_f(0) = (((r_in + x_lo*r_diff)*sin(2.*PI*y_hi))
				                                  -((r_in + x_lo*r_diff)*sin(2.*PI*y_lo)))/a_dy;
			if (a_s == 1 && a_d == 0) a_N_ave_f(0) = -(((r_in + x_lo*r_diff)*cos(2.*PI*y_hi))
				                                   -((r_in + x_lo*r_diff)*cos(2.*PI*y_lo)))/a_dy;
			if (a_s == 0 && a_d == 1) a_N_ave_f(0) = -(((r_in + x_hi*r_diff)*sin(2.*PI*y_lo))
				                                   -((r_in + x_lo*r_diff)*sin(2.*PI*y_lo)))/a_dx;
			if (a_s == 1 && a_d == 1) a_N_ave_f(0) = (((r_in + x_hi*r_diff)*cos(2.*PI*y_lo))
				                                  -((r_in + x_lo*r_diff)*cos(2.*PI*y_lo)))/a_dx;
		}
#endif

#if DIM == 3
		if (grid_type_global == 0) {
			if (a_s == 0 && a_d == 0) a_N_ave_f(0) = 1.0;
			if (a_s == 1 && a_d == 0) a_N_ave_f(0) = 0.0;
			if (a_s == 2 && a_d == 0) a_N_ave_f(0) = 0.0;
			if (a_s == 0 && a_d == 1) a_N_ave_f(0) = 0.0;
			if (a_s == 1 && a_d == 1) a_N_ave_f(0) = 1.0;
			if (a_s == 2 && a_d == 1) a_N_ave_f(0) = 0.0;
			if (a_s == 0 && a_d == 2) a_N_ave_f(0) = 0.0;
			if (a_s == 1 && a_d == 2) a_N_ave_f(0) = 0.0;
			if (a_s == 2 && a_d == 2) a_N_ave_f(0) = 1.0;
		}

		if (grid_type_global == 2) {
			double eta1_hi = (a_pt[0] + 1.0)*a_dx;
			double eta1_lo = (a_pt[0] - 0.0)*a_dx;
			double eta2_hi = (a_pt[1] + 1.0)*a_dy;
			double eta2_lo = (a_pt[1] - 0.0)*a_dy;
			double eta3_hi = (a_pt[2] + 1.0)*a_dz;
			double eta3_lo = (a_pt[2] - 0.0)*a_dz;
			double R_t = (r_out - r_in)/(exp(C_rad) - 1.0);
			double r_lo = r_in + R_t*(exp(C_rad*eta1_lo) - 1.0);

			auto NN_21_1 = [](double a_R_t, double eta1,double eta2,double eta3){
				double r = r_in + a_R_t*(exp(C_rad*eta1) - 1.0);
				return PI*r*r*sin(2.0*PI*eta2)*eta3/2.0;
			};

			auto NN_21_2 = [](double a_R_t, double eta1,double eta2,double eta3){
				double r = r_in + a_R_t*(exp(C_rad*eta1) - 1.0);
				return -PI*r*r*cos(2.0*PI*eta2)*eta3/2.0;
			};

			auto NN_31_1 = [](double a_R_t, double eta1,double eta2,double eta3){
				double r = r_in + a_R_t*(exp(C_rad*eta1) - 1.0);
				return -0.5*r*r*sin(2.0*PI*eta2)*sin(PI*eta3)*cos(PI*eta3);
			};

			auto NN_31_2 = [](double a_R_t, double eta1,double eta2,double eta3){
				double r = r_in + a_R_t*(exp(C_rad*eta1) - 1.0);
				return 0.5*r*r*cos(2.0*PI*eta2)*sin(PI*eta3)*cos(PI*eta3);
			};

			auto NN_31_3 = [](double a_R_t, double eta1,double eta2,double eta3){
				double r = r_in + a_R_t*(exp(C_rad*eta1) - 1.0);
				return PI*r*r*sin(PI*eta3)*sin(PI*eta3)*eta2;
			};
			if (a_s == 0 && a_d == 0) a_N_ave_f(0) =(+(NN_21_1(R_t,eta1_lo,eta2_hi,eta3_hi)-NN_21_1(R_t,eta1_lo,eta2_hi,eta3_lo))
													 -(NN_31_1(R_t,eta1_lo,eta2_lo,eta3_hi)-NN_31_1(R_t,eta1_lo,eta2_hi,eta3_hi))
													 +(NN_21_1(R_t,eta1_lo,eta2_lo,eta3_lo)-NN_21_1(R_t,eta1_lo,eta2_lo,eta3_hi))
													 -(NN_31_1(R_t,eta1_lo,eta2_hi,eta3_lo)-NN_31_1(R_t,eta1_lo,eta2_lo,eta3_lo)))/a_dy/a_dz;
			if (a_s == 1 && a_d == 0) a_N_ave_f(0) =(+(NN_21_2(R_t,eta1_lo,eta2_hi,eta3_hi)-NN_21_2(R_t,eta1_lo,eta2_hi,eta3_lo))
													 -(NN_31_2(R_t,eta1_lo,eta2_lo,eta3_hi)-NN_31_2(R_t,eta1_lo,eta2_hi,eta3_hi))
													 +(NN_21_2(R_t,eta1_lo,eta2_lo,eta3_lo)-NN_21_2(R_t,eta1_lo,eta2_lo,eta3_hi))
													 -(NN_31_2(R_t,eta1_lo,eta2_hi,eta3_lo)-NN_31_2(R_t,eta1_lo,eta2_lo,eta3_lo)))/a_dy/a_dz;
			if (a_s == 2 && a_d == 0) a_N_ave_f(0) =(-(NN_31_3(R_t,eta1_lo,eta2_lo,eta3_hi)-NN_31_3(R_t,eta1_lo,eta2_hi,eta3_hi))
													 -(NN_31_3(R_t,eta1_lo,eta2_hi,eta3_lo)-NN_31_3(R_t,eta1_lo,eta2_lo,eta3_lo)))/a_dy/a_dz;
			if (a_s == 0 && a_d == 1) a_N_ave_f(0) =(+(NN_21_1(R_t,eta1_hi,eta2_lo,eta3_lo)-NN_21_1(R_t,eta1_hi,eta2_lo,eta3_hi))
													 +(NN_21_1(R_t,eta1_lo,eta2_lo,eta3_hi)-NN_21_1(R_t,eta1_lo,eta2_lo,eta3_lo)))/a_dx/a_dz;
			if (a_s == 1 && a_d == 1) a_N_ave_f(0) =(+(NN_21_2(R_t,eta1_hi,eta2_lo,eta3_lo)-NN_21_2(R_t,eta1_hi,eta2_lo,eta3_hi))
													 +(NN_21_2(R_t,eta1_lo,eta2_lo,eta3_hi)-NN_21_2(R_t,eta1_lo,eta2_lo,eta3_lo)))/a_dx/a_dz;
			if (a_s == 2 && a_d == 1) a_N_ave_f(0) = 0.0;
			if (a_s == 0 && a_d == 2) a_N_ave_f(0) =(-(NN_31_1(R_t,eta1_hi,eta2_hi,eta3_lo)-NN_31_1(R_t,eta1_hi,eta2_lo,eta3_lo))
													 -(NN_31_1(R_t,eta1_lo,eta2_lo,eta3_lo)-NN_31_1(R_t,eta1_lo,eta2_hi,eta3_lo)))/a_dx/a_dy;
			if (a_s == 1 && a_d == 2) a_N_ave_f(0) =(-(NN_31_2(R_t,eta1_hi,eta2_hi,eta3_lo)-NN_31_2(R_t,eta1_hi,eta2_lo,eta3_lo))
													 -(NN_31_2(R_t,eta1_lo,eta2_lo,eta3_lo)-NN_31_2(R_t,eta1_lo,eta2_hi,eta3_lo)))/a_dx/a_dy;
			if (a_s == 2 && a_d == 2) a_N_ave_f(0) =(-(NN_31_3(R_t,eta1_hi,eta2_hi,eta3_lo)-NN_31_3(R_t,eta1_hi,eta2_lo,eta3_lo))
													 -(NN_31_3(R_t,eta1_lo,eta2_lo,eta3_lo)-NN_31_3(R_t,eta1_lo,eta2_hi,eta3_lo)))/a_dx/a_dy;


		}
#endif
	}
	PROTO_KERNEL_END(N_ave_f_calcF, N_ave_f_calc)





	PROTO_KERNEL_START
	void N_ave_f_calc2F( const Point& a_pt,
	                    Var<double,DIM*DIM>& a_N_ave_f,
	                    const double a_dx,
	                    const double a_dy,
	                    const double a_dz)
	{
#if DIM == 2
		int x_loc = a_pt[0];
		int y_loc = a_pt[1];
		double x_hi = (x_loc + 1.0)*a_dx;
		double x_lo = (x_loc - 0.0)*a_dx;
		double y_hi = (y_loc + 1.0)*a_dy;
		double y_lo = (y_loc - 0.0)*a_dy;

		if (grid_type_global == 0) {
			a_N_ave_f(0) = 1.0;
			a_N_ave_f(1) = 0.0;
			a_N_ave_f(2) = 0.0;
			a_N_ave_f(3) = 1.0;
		}

		if (grid_type_global == 1) {
			double C1 = C1_fix;
			double C2 = C1;
			a_N_ave_f(0) = ((y_hi+C2*sin(2.*PI*x_lo)*sin(2.*PI*y_hi))
				                                  -(y_lo+C2*sin(2.*PI*x_lo)*sin(2.*PI*y_lo)))/a_dy;
			a_N_ave_f(1) = -((x_lo+C1*sin(2.*PI*x_lo)*sin(2.*PI*y_hi))
				                                   -(x_lo+C1*sin(2.*PI*x_lo)*sin(2.*PI*y_lo)))/a_dy;
			a_N_ave_f(2) = -((y_lo+C2*sin(2.*PI*x_hi)*sin(2.*PI*y_lo))
				                                   -(y_lo+C2*sin(2.*PI*x_lo)*sin(2.*PI*y_lo)))/a_dx;
			a_N_ave_f(3) = ((x_hi+C1*sin(2.*PI*x_hi)*sin(2.*PI*y_lo))
				                                  -(x_lo+C1*sin(2.*PI*x_lo)*sin(2.*PI*y_lo)))/a_dx;
		}

		if (grid_type_global == 2) {
			double r_diff = r_out-r_in;
			a_N_ave_f(0) = (((r_in + x_lo*r_diff)*sin(2.*PI*y_hi))
				                                  -((r_in + x_lo*r_diff)*sin(2.*PI*y_lo)))/a_dy;
			a_N_ave_f(1) = -(((r_in + x_lo*r_diff)*cos(2.*PI*y_hi))
				                                   -((r_in + x_lo*r_diff)*cos(2.*PI*y_lo)))/a_dy;
			a_N_ave_f(2) = -(((r_in + x_hi*r_diff)*sin(2.*PI*y_lo))
				                                   -((r_in + x_lo*r_diff)*sin(2.*PI*y_lo)))/a_dx;
			a_N_ave_f(3) = (((r_in + x_hi*r_diff)*cos(2.*PI*y_lo))
				                                  -((r_in + x_lo*r_diff)*cos(2.*PI*y_lo)))/a_dx;
		}
#endif

#if DIM == 3
		if (grid_type_global == 0) {
			a_N_ave_f(0) = 1.0;
			a_N_ave_f(1) = 0.0;
			a_N_ave_f(2) = 0.0;
			a_N_ave_f(3) = 0.0;
			a_N_ave_f(4) = 1.0;
			a_N_ave_f(5) = 0.0;
			a_N_ave_f(6) = 0.0;
			a_N_ave_f(7) = 0.0;
			a_N_ave_f(8) = 1.0;
		}

		if (grid_type_global == 2) {
			double eta1_hi = (a_pt[0] + 1.0)*a_dx;
			double eta1_lo = (a_pt[0] - 0.0)*a_dx;
			double eta2_hi = (a_pt[1] + 1.0)*a_dy;
			double eta2_lo = (a_pt[1] - 0.0)*a_dy;
			double eta3_hi = (a_pt[2] + 1.0)*a_dz;
			double eta3_lo = (a_pt[2] - 0.0)*a_dz;
			double R_t = (r_out - r_in)/(exp(C_rad) - 1.0);
			double r_lo = r_in + R_t*(exp(C_rad*eta1_lo) - 1.0);

			auto NN_21_1 = [](double a_R_t, double eta1,double eta2,double eta3){
				double r = r_in + a_R_t*(exp(C_rad*eta1) - 1.0);
				return PI*r*r*sin(2.0*PI*eta2)*eta3/2.0;
			};

			auto NN_21_2 = [](double a_R_t, double eta1,double eta2,double eta3){
				double r = r_in + a_R_t*(exp(C_rad*eta1) - 1.0);
				return -PI*r*r*cos(2.0*PI*eta2)*eta3/2.0;
			};

			auto NN_31_1 = [](double a_R_t, double eta1,double eta2,double eta3){
				double r = r_in + a_R_t*(exp(C_rad*eta1) - 1.0);
				return -0.5*r*r*sin(2.0*PI*eta2)*sin(PI*eta3)*cos(PI*eta3);
			};

			auto NN_31_2 = [](double a_R_t, double eta1,double eta2,double eta3){
				double r = r_in + a_R_t*(exp(C_rad*eta1) - 1.0);
				return 0.5*r*r*cos(2.0*PI*eta2)*sin(PI*eta3)*cos(PI*eta3);
			};

			auto NN_31_3 = [](double a_R_t, double eta1,double eta2,double eta3){
				double r = r_in + a_R_t*(exp(C_rad*eta1) - 1.0);
				return PI*r*r*sin(PI*eta3)*sin(PI*eta3)*eta2;
			};
			a_N_ave_f(0) =(+(NN_21_1(R_t,eta1_lo,eta2_hi,eta3_hi)-NN_21_1(R_t,eta1_lo,eta2_hi,eta3_lo))
													 -(NN_31_1(R_t,eta1_lo,eta2_lo,eta3_hi)-NN_31_1(R_t,eta1_lo,eta2_hi,eta3_hi))
													 +(NN_21_1(R_t,eta1_lo,eta2_lo,eta3_lo)-NN_21_1(R_t,eta1_lo,eta2_lo,eta3_hi))
													 -(NN_31_1(R_t,eta1_lo,eta2_hi,eta3_lo)-NN_31_1(R_t,eta1_lo,eta2_lo,eta3_lo)))/a_dy/a_dz;
			a_N_ave_f(1) =(+(NN_21_2(R_t,eta1_lo,eta2_hi,eta3_hi)-NN_21_2(R_t,eta1_lo,eta2_hi,eta3_lo))
													 -(NN_31_2(R_t,eta1_lo,eta2_lo,eta3_hi)-NN_31_2(R_t,eta1_lo,eta2_hi,eta3_hi))
													 +(NN_21_2(R_t,eta1_lo,eta2_lo,eta3_lo)-NN_21_2(R_t,eta1_lo,eta2_lo,eta3_hi))
													 -(NN_31_2(R_t,eta1_lo,eta2_hi,eta3_lo)-NN_31_2(R_t,eta1_lo,eta2_lo,eta3_lo)))/a_dy/a_dz;
			a_N_ave_f(2) =(-(NN_31_3(R_t,eta1_lo,eta2_lo,eta3_hi)-NN_31_3(R_t,eta1_lo,eta2_hi,eta3_hi))
													 -(NN_31_3(R_t,eta1_lo,eta2_hi,eta3_lo)-NN_31_3(R_t,eta1_lo,eta2_lo,eta3_lo)))/a_dy/a_dz;
			a_N_ave_f(3) =(+(NN_21_1(R_t,eta1_hi,eta2_lo,eta3_lo)-NN_21_1(R_t,eta1_hi,eta2_lo,eta3_hi))
													 +(NN_21_1(R_t,eta1_lo,eta2_lo,eta3_hi)-NN_21_1(R_t,eta1_lo,eta2_lo,eta3_lo)))/a_dx/a_dz;
			a_N_ave_f(4) =(+(NN_21_2(R_t,eta1_hi,eta2_lo,eta3_lo)-NN_21_2(R_t,eta1_hi,eta2_lo,eta3_hi))
													 +(NN_21_2(R_t,eta1_lo,eta2_lo,eta3_hi)-NN_21_2(R_t,eta1_lo,eta2_lo,eta3_lo)))/a_dx/a_dz;
			a_N_ave_f(5) = 0.0;
			a_N_ave_f(6) =(-(NN_31_1(R_t,eta1_hi,eta2_hi,eta3_lo)-NN_31_1(R_t,eta1_hi,eta2_lo,eta3_lo))
													 -(NN_31_1(R_t,eta1_lo,eta2_lo,eta3_lo)-NN_31_1(R_t,eta1_lo,eta2_hi,eta3_lo)))/a_dx/a_dy;
			a_N_ave_f(7) =(-(NN_31_2(R_t,eta1_hi,eta2_hi,eta3_lo)-NN_31_2(R_t,eta1_hi,eta2_lo,eta3_lo))
													 -(NN_31_2(R_t,eta1_lo,eta2_lo,eta3_lo)-NN_31_2(R_t,eta1_lo,eta2_hi,eta3_lo)))/a_dx/a_dy;
			a_N_ave_f(8) =(-(NN_31_3(R_t,eta1_hi,eta2_hi,eta3_lo)-NN_31_3(R_t,eta1_hi,eta2_lo,eta3_lo))
													 -(NN_31_3(R_t,eta1_lo,eta2_lo,eta3_lo)-NN_31_3(R_t,eta1_lo,eta2_hi,eta3_lo)))/a_dx/a_dy;


		}
#endif
	}
	PROTO_KERNEL_END(N_ave_f_calc2F, N_ave_f_calc2)







	PROTO_KERNEL_START
	void d_U_calcF(  State& a_d_U,
	                 const State& a_UJ_ahead,
	                 const State& a_UJ_behind,
	                 const Var<double,1>& a_J_ahead,
	                 const Var<double,1>& a_J_behind)
	{
		for (int i=0; i< NUMCOMPS; i++) {
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
		for (int i=0; i< NUMCOMPS; i++) {
			a_U_demapped(i) = (a_a_U(i)-(a_dot_pro_sum2(i)/12.0))/a_Jacobian_ave(0);
		}

	}
	PROTO_KERNEL_END(a_U_demapped_calcF, a_U_demapped_calc)


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


	void N_ave_f_calc_func(BoxData<double,1>& a_N_ave_f,
	                       const int a_s,
	                       const int a_d,
	                       const double a_dx,
	                       const double a_dy,
	                       const double a_dz)
	{
		forallInPlace_p(N_ave_f_calc, a_N_ave_f, a_s, a_d, a_dx, a_dy, a_dz);
	}

	void N_ave_f_calc_func(BoxData<double,DIM*DIM>& a_N_ave_f,
	                       const double a_dx,
	                       const double a_dy,
	                       const double a_dz)
	{
		forallInPlace_p(N_ave_f_calc2, a_N_ave_f, a_dx, a_dy, a_dz);
	}


	void Jacobian_Ave_calc(BoxData<double,1>& a_Jacobian_ave,
	                       const double a_dx,
	                       const double a_dy,
	                       const double a_dz,
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
		double dxd[3] = {a_dx, a_dy, a_dz};
		a_Jacobian_ave.setVal(0.0);

		for (int d = 0; d < DIM; d++) {
			Scalar X_f_mapped(a_dbx0);
			X_f_mapped.setVal(0.0);
			for (int s = 0; s < DIM; s++) {
				forallInPlace_p(N_ave_f_calc, N_s_d_ave_f, s, d, a_dx, a_dy, a_dz);
				forallInPlace_p(X_ave_f_calc, X_ave_f, s, d, a_dx, a_dy, a_dz);

				// if (s==0 && d==0){
				// std::string filename="N_s_d_ave_f";
				// WriteBoxData(filename.c_str(),N_s_d_ave_f,a_dx);
				// filename="X_ave_f";
				// WriteBoxData(filename.c_str(),X_ave_f,a_dx);
				// }

				Scalar dot_pro_sum(a_dbx0);
				dot_pro_sum.setVal(0.0);
				for (int s_temp = 0; s_temp < DIM; s_temp++) {
					if (s_temp != d) {
						Scalar d_perp_N_s = m_derivative[s_temp](N_s_d_ave_f);
						Scalar d_perp_X = m_derivative[s_temp](X_ave_f);
						Scalar dot_pro = forall<double,1>(dot_pro_calcX,d_perp_N_s,d_perp_X);
						dot_pro_sum += dot_pro;
					}
				}
				Scalar X_f_mapped1D = forall<double,1>(X_f_mapped1D_calc,X_ave_f,N_s_d_ave_f,dot_pro_sum);
				X_f_mapped += X_f_mapped1D;
			}
			Scalar Rhs_d = m_divergence[d](X_f_mapped);
			Rhs_d *= 1./dxd[d];
			a_Jacobian_ave += Rhs_d;
		}
		// std::string filename="a_Jacobian_ave";
		// WriteBoxData(filename.c_str(),a_Jacobian_ave,a_dx);

	}



	void JU_to_U_calc(BoxData<double,NUMCOMPS>& a_U_demapped,
	                  const BoxData<double,NUMCOMPS>& a_a_U,
	                  BoxData<double,1>& a_Jacobian_ave,
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
		for (int d = 0; d < DIM; d++) {

			Vector UJ_ahead = m_ahead_shift[d](a_a_U);
			Vector UJ_behind = m_behind_shift[d](a_a_U);
			Scalar J_ahead = alias(a_Jacobian_ave,Point::Basis(d)*(-1));
			Scalar J_behind = alias(a_Jacobian_ave,Point::Basis(d)*(1));

			Vector d_U = forall<double,NUMCOMPS>(d_U_calc,UJ_ahead,UJ_behind,J_ahead,J_behind);
			Scalar d_J = m_derivative[d](a_Jacobian_ave);
			Vector dot_pro = forall<double,NUMCOMPS>(dot_pro_calcF,d_J,d_U);
			dot_pro_sum2 += dot_pro;
		}

		a_U_demapped = forall<double,NUMCOMPS>(a_U_demapped_calc,a_Jacobian_ave,a_a_U,dot_pro_sum2);
	}


	PROTO_KERNEL_START
	void a_U_demapped_2nd_order_calcF(  State& a_U_demapped,
	                                    const Var<double,1>& a_Jacobian_ave,
	                                    const State& a_a_U)
	{
		for (int i=0; i< NUMCOMPS; i++) {
			a_U_demapped(i) = a_a_U(i)/a_Jacobian_ave(0);
		}

	}
	PROTO_KERNEL_END(a_U_demapped_2nd_order_calcF, a_U_demapped_2nd_order_calc)


	void JU_to_U_2ndOrdercalc(BoxData<double,NUMCOMPS>& a_U_demapped,
	                          const BoxData<double,NUMCOMPS>& a_a_U,
	                          const BoxData<double,1>& a_Jacobian_ave,
	                          const Box& a_dbx0)
	{
		a_U_demapped = forall<double,NUMCOMPS>(a_U_demapped_2nd_order_calc,a_Jacobian_ave,a_a_U);
	}


	void JU_to_W_calc(BoxData<double,NUMCOMPS>& a_W,
	                  const BoxData<double,NUMCOMPS>& a_JU,
	                  const double a_dx,
	                  const double a_dy,
	                  const double a_dz,
	                  const double a_gamma)
	{
		Box dbx0 = a_JU.box();
		static Stencil<double> m_laplacian;
		static Stencil<double> m_deconvolve;
		static Stencil<double> m_copy;
		static bool initialized = false;
		if(!initialized)
		{
			m_laplacian = Stencil<double>::Laplacian();
			m_deconvolve = (-1.0/24.0)*m_laplacian + (1.0)*Shift(Point::Zeros());
			m_copy = 1.0*Shift(Point::Zeros());
			initialized =  true;
		}


		using namespace std;
		a_W.setVal(0.0);
		double gamma = a_gamma;

		Scalar Jacobian_ave(dbx0);
		MHD_Mapping::Jacobian_Ave_calc(Jacobian_ave,a_dx, a_dy, a_dz,dbx0);

		Vector a_U(dbx0);

		MHD_Mapping::JU_to_U_calc(a_U, a_JU, Jacobian_ave, dbx0);

		Vector W_bar(dbx0);
		MHDOp::consToPrimcalc(W_bar,a_U,gamma);
		Vector U = m_deconvolve(a_U);
		Vector W(dbx0);
		MHDOp::consToPrimcalc(W,U,gamma);
		a_W = m_laplacian(W_bar,1.0/24.0);
		a_W += W;

	}




	void JU_to_W_bar_calc(BoxData<double,NUMCOMPS>& a_W_bar,
	                      const BoxData<double,NUMCOMPS>& a_JU,
	                      const double a_dx,
	                      const double a_dy,
	                      const double a_dz,
	                      const double a_gamma)
	{
		Box dbx1 = a_JU.box();
		Box dbx0 = dbx1.grow(NGHOST);
		a_W_bar.setVal(0.0);
		double gamma = a_gamma;
		Scalar Jacobian_ave(dbx0);
		MHD_Mapping::Jacobian_Ave_calc(Jacobian_ave,a_dx, a_dy, a_dz,dbx0);
		Vector a_U(dbx0);
		MHD_Mapping::JU_to_U_calc(a_U, a_JU, Jacobian_ave, dbx0);
		MHDOp::consToPrimcalc(a_W_bar,a_U,gamma);
	}



	PROTO_KERNEL_START
	void out_data_joinF(Var<double,NUMCOMPS+DIM>& a_out_data,
	                    const Var<double,DIM>& a_phys_coords,
	                    const Var<double,NUMCOMPS>& a_W)
	{
		for (int i=0; i< DIM; i++) {
			a_out_data(i) = a_phys_coords(i);
		}
		for (int i=DIM; i< NUMCOMPS+DIM; i++) {
			a_out_data(i) = a_W(i-DIM);
		}

	}
	PROTO_KERNEL_END(out_data_joinF, out_data_join)


	void out_data_calc(BoxData<double,NUMCOMPS+DIM>& a_out_data,
	                   const BoxData<double,DIM>& a_phys_coords,
	                   const BoxData<double,NUMCOMPS>& a_W)
	{
		a_out_data = forall<double,NUMCOMPS+DIM>(out_data_join, a_phys_coords, a_W);
	}

	PROTO_KERNEL_START
	void X_ave_calcF( const Point& a_pt,
	                    Var<double,DIM>& a_X_ave,
	                    const double a_dx,
	                    const double a_dy,
	                    const double a_dz)
	{

#if DIM == 2
		double eta1_hi = (a_pt[0] + 1.0)*a_dx;
		double eta1_lo = (a_pt[0] - 0.0)*a_dx;
		double eta2_hi = (a_pt[1] + 1.0)*a_dy;
		double eta2_lo = (a_pt[1] - 0.0)*a_dy;


		if (grid_type_global == 2) {
			double r_diff = r_out-r_in;
			a_X_ave(0) = (r_in*eta1_hi + r_diff*eta1_hi*eta1_hi/2.0 - r_in*eta1_lo - r_diff*eta1_lo*eta1_lo/2.0)*(sin(2.0*PI*eta2_hi)-sin(2.0*PI*eta2_lo))/a_dx/a_dy/2.0/PI;
			a_X_ave(1) = (r_in*eta1_hi + r_diff*eta1_hi*eta1_hi/2.0 - r_in*eta1_lo - r_diff*eta1_lo*eta1_lo/2.0)*(cos(2.0*PI*eta2_lo)-cos(2.0*PI*eta2_hi))/a_dx/a_dy/2.0/PI;

		}
#endif

#if DIM == 3

		if (grid_type_global == 2) {
			double eta1_hi = (a_pt[0] + 1.0)*a_dx;
			double eta1_lo = (a_pt[0] - 0.0)*a_dx;
			double eta2_hi = (a_pt[1] + 1.0)*a_dy;
			double eta2_lo = (a_pt[1] - 0.0)*a_dy;
			double eta3_hi = (a_pt[2] + 1.0)*a_dz;
			double eta3_lo = (a_pt[2] - 0.0)*a_dz;

			double R_t = (r_out - r_in)/(exp(C_rad) - 1.0);
			double r_lo = r_in + R_t*(exp(C_rad*eta1_lo) - 1.0);
			double coseta3_hi = cos(PI*eta3_hi);
			double coseta3_lo = cos(PI*eta3_lo);
			double coseta2_hi = cos(2.0*PI*eta2_hi);
			double coseta2_lo = cos(2.0*PI*eta2_lo);
			double sineta3_hi = sin(PI*eta3_hi);
			double sineta3_lo = sin(PI*eta3_lo);
			double sineta2_hi = sin(2.0*PI*eta2_hi);
			double sineta2_lo = sin(2.0*PI*eta2_lo);

			double r_int_lo = r_in*eta1_lo + R_t*((exp(C_rad*eta1_lo)/C_rad) - eta1_lo);
			double r_int_hi = r_in*eta1_hi + R_t*((exp(C_rad*eta1_hi)/C_rad) - eta1_hi);

			// To be finished


		}

#endif

	}
	PROTO_KERNEL_END(X_ave_calcF, X_ave_calc)

	PROTO_KERNEL_START
	void Cart_to_Sph_velocityF(Var<double,NUMCOMPS>& a_W_Sph,
	                    const Var<double,DIM>& a_X_ave,
	                    const Var<double,NUMCOMPS>& a_W_Cart)
	{
#if DIM == 2
		//double r = sqrt(a_X_ave(0)*a_X_ave(0) + a_X_ave(1)*a_X_ave(1));
		double theta = atan2(a_X_ave(1),a_X_ave(0));
		double u_r =  a_W_Cart(1)*cos(theta) + a_W_Cart(2)*sin(theta);
		double u_t = -a_W_Cart(1)*sin(theta) + a_W_Cart(2)*cos(theta);

		a_W_Sph(0) = a_W_Cart(0);
		a_W_Sph(1) = u_r;
		a_W_Sph(2) = u_t;
		//a_W_Sph(2) = 0.0;
		a_W_Sph(3) = a_W_Cart(3);
		a_W_Sph(4) = a_W_Cart(4);
		a_W_Sph(5) = a_W_Cart(5);
#endif
	}
	PROTO_KERNEL_END(Cart_to_Sph_velocityF, Cart_to_Sph_velocity)

	void W_Cart_to_W_Sph(BoxData<double,NUMCOMPS>& a_W_Sph,
	                      const BoxData<double,NUMCOMPS>& a_W_Cart,
	                      const double a_dx,
	                      const double a_dy,
	                      const double a_dz)
	{
		Box dbx0 = a_W_Cart.box();
		BoxData<double,DIM,1,1> X_ave(dbx0);
		forallInPlace_p(X_ave_calc, X_ave, a_dx, a_dy, a_dz);
		a_W_Sph = forall<double,NUMCOMPS>(Cart_to_Sph_velocity, X_ave, a_W_Cart);
	}


	PROTO_KERNEL_START
	void Cart_to_Sph_velocity_edgeF(Var<double,NUMCOMPS>& a_W_Sph,
									const Var<double,1>& a_X_ave_f,
									const Var<double,1>& a_Y_ave_f,
									const Var<double,1>& a_Z_ave_f,
									const Var<double,NUMCOMPS>& a_W_Cart)
	{
#if DIM == 2
		//double r = sqrt(a_X_ave(0)*a_X_ave(0) + a_X_ave(1)*a_X_ave(1));
		double theta = atan2(a_Y_ave_f(0),a_X_ave_f(0));
		double u_r =  a_W_Cart(1)*cos(theta) + a_W_Cart(2)*sin(theta);
		double u_t = -a_W_Cart(1)*sin(theta) + a_W_Cart(2)*cos(theta);

		a_W_Sph(0) = a_W_Cart(0);
		a_W_Sph(1) = u_r;
		a_W_Sph(2) = u_t;
		a_W_Sph(3) = a_W_Cart(3);
		a_W_Sph(4) = a_W_Cart(4);
		a_W_Sph(5) = a_W_Cart(5);
#endif
	}
	PROTO_KERNEL_END(Cart_to_Sph_velocity_edgeF, Cart_to_Sph_velocity_edge)


	void W_Cart_to_W_Sph_edge(BoxData<double,NUMCOMPS>& a_W_Sph,
	                      const BoxData<double,NUMCOMPS>& a_W_Cart,
						  int a_d,
	                      const double a_dx,
	                      const double a_dy,
	                      const double a_dz)
	{
		Box dbx0 = a_W_Cart.box();
		Scalar X_ave_f(dbx0), Y_ave_f(dbx0), Z_ave_f(dbx0);
		int s = 0;
		forallInPlace_p(X_ave_f_calc, X_ave_f, s, a_d, a_dx, a_dy, a_dz);
		s = 1;
		forallInPlace_p(X_ave_f_calc, Y_ave_f, s, a_d, a_dx, a_dy, a_dz);
		s = 2;
		forallInPlace_p(X_ave_f_calc, Z_ave_f, s, a_d, a_dx, a_dy, a_dz);
		a_W_Sph = forall<double,NUMCOMPS>(Cart_to_Sph_velocity_edge, X_ave_f, Y_ave_f, Z_ave_f, a_W_Cart);
	}


	PROTO_KERNEL_START
	void Sph_to_Cart_velocityF(Var<double,NUMCOMPS>& a_W_Cart,
	                    const Var<double,1>& a_X_ave_f,
	                    const Var<double,1>& a_Y_ave_f,
	                    const Var<double,1>& a_Z_ave_f,
	                    const Var<double,NUMCOMPS>& a_W_Sph)
	{
#if DIM == 2
		double theta = atan2(a_Y_ave_f(0),a_X_ave_f(0));
		double u_x = a_W_Sph(1)*cos(theta) - a_W_Sph(2)*sin(theta);
		double u_y = a_W_Sph(1)*sin(theta) + a_W_Sph(2)*cos(theta);

		a_W_Cart(0) = a_W_Sph(0);
		a_W_Cart(1) = u_x;
		a_W_Cart(2) = u_y;
		a_W_Cart(3) = a_W_Sph(3);
		a_W_Cart(4) = a_W_Sph(4);
		a_W_Cart(5) = a_W_Sph(5);
#endif
	}
	PROTO_KERNEL_END(Sph_to_Cart_velocityF, Sph_to_Cart_velocity)




	void W_Sph_to_W_Cart(BoxData<double,NUMCOMPS>& a_W_Cart,
	                      const BoxData<double,NUMCOMPS>& a_W_Sph,
						  int a_d,
	                      const double a_dx,
	                      const double a_dy,
	                      const double a_dz)
	{
		Box dbx0 = a_W_Sph.box();
		Scalar X_ave_f(dbx0), Y_ave_f(dbx0), Z_ave_f(dbx0);
		int s = 0;
		forallInPlace_p(X_ave_f_calc, X_ave_f, s, a_d, a_dx, a_dy, a_dz);
		s = 1;
		forallInPlace_p(X_ave_f_calc, Y_ave_f, s, a_d, a_dx, a_dy, a_dz);
		s = 2;
		forallInPlace_p(X_ave_f_calc, Z_ave_f, s, a_d, a_dx, a_dy, a_dz);
		a_W_Cart = forall<double,NUMCOMPS>(Sph_to_Cart_velocity, X_ave_f, Y_ave_f, Z_ave_f, a_W_Sph);

	}


}
