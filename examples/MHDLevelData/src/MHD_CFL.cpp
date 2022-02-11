#include "Proto.H"
#include "MHD_Mapping.H"
#include "MHDOp.H"
#include "MHD_Initialize.H"
#include "MHD_Output_Writer.H"
#include "Proto_WriteBoxData.H"
#include <algorithm>    // std::min

extern int init_condition_type;
extern int grid_type_global;
extern double r_in;
extern double r_out;
extern bool initialize_in_spherical_coords;

typedef BoxData<double,1> Scalar;
typedef BoxData<double,NUMCOMPS> Vector;

namespace MHD_CFL {

    PROTO_KERNEL_START
	void lambdacalcF(Var<double,1>& a_lambda,
	                 const State& a_W_ave,
	                 int a_d,
	                 double a_gamma)
	{
		double gamma = a_gamma;
		double rho=0., u=0., v=0., w=0., p=0., Bx=0., By=0., Bz=0., ce, af, B_mag, u_mag, Bdir, udir;
#if DIM == 1
		rho = a_W_ave(0);
		u   = a_W_ave(1);
		p   = a_W_ave(2);
		Bx  = a_W_ave(3);
#endif
#if DIM == 2
		rho = a_W_ave(0);
		u   = a_W_ave(1);
		v   = a_W_ave(2);
		p   = a_W_ave(3);
		Bx  = a_W_ave(4);
		By  = a_W_ave(5);
#endif
#if DIM == 3
		rho = a_W_ave(0);
		u   = a_W_ave(1);
		v   = a_W_ave(2);
		w   = a_W_ave(3);
		p   = a_W_ave(4);
		Bx  = a_W_ave(5);
		By  = a_W_ave(6);
		Bz  = a_W_ave(7);
#endif
		if (a_d == 0) {
			Bdir = Bx;
			udir = u;
		};
		if (a_d == 1) {
			Bdir = By;
			udir = v;
		};
		if (a_d == 2) {
			Bdir = Bz;
			udir = w;
		};

		if (p < 0.0) p = 0.0;
		ce = sqrt(gamma*p/rho);
		B_mag = sqrt(Bx*Bx+By*By+Bz*Bz);
		af = sqrt(ce*ce + B_mag*B_mag/4.0/PI/rho);
		u_mag = sqrt(u*u+v*v+w*w);
		// a_lambda(0) = af + u_mag;
		a_lambda(0) = af + abs(udir);

	}
	PROTO_KERNEL_END(lambdacalcF, lambdacalc)


    PROTO_KERNEL_START
	void dt_dcalcF(Var<double,1>& a_dt_d,
	               Var<double,1>& a_Lambda_f,
                   V& a_x,
                   V& a_x_ahead)
	{
        double dx_d;
        #if DIM == 2
            dx_d = sqrt((a_x_ahead(0)-a_x(0))*(a_x_ahead(0)-a_x(0)) + (a_x_ahead(1)-a_x(1))*(a_x_ahead(1)-a_x(1)));
        #endif
        #if DIM == 3
            dx_d = sqrt((a_x_ahead(0)-a_x(0))*(a_x_ahead(0)-a_x(0)) + (a_x_ahead(1)-a_x(1))*(a_x_ahead(1)-a_x(1)) + (a_x_ahead(2)-a_x(2))*(a_x_ahead(2)-a_x(2)));
        #endif
		a_dt_d(0) = dx_d/a_Lambda_f(0);
	}
	PROTO_KERNEL_END(dt_dcalcF, dt_dcalc)


	void Min_dt_calc_func(double& a_dt,
	                    const BoxData<double,NUMCOMPS>& a_W_ave,
                        const double a_dx,
                        const double a_dy,
                        const double a_dz,
	                    const double a_gamma)
	{   
        double dt[DIM];
        for (int dir = 0; dir < DIM; dir++)
		{
		    Scalar Lambda_f = forall<double>(lambdacalc, a_W_ave, dir, a_gamma);
            Box dbx0 = a_W_ave.box();
            BoxData<double,DIM> eta(dbx0);
            MHD_Mapping::etaFace_calc(eta, dbx0, a_dx, a_dy, a_dz, dir);
            BoxData<double,DIM> x(dbx0);		
            MHD_Mapping::eta_to_x_calc(x,eta);
            BoxData<double,DIM> x_ahead = alias(x,Point::Basis(dir)*(-1));
            Scalar dt_d = forall<double>(dt_dcalc, Lambda_f, x, x_ahead);
			// std::string filename="dt_d";
			// if (dir == 0) WriteBoxData(filename.c_str(),dt_d,a_dx);
            dt[dir] = dt_d.min();
        }
        #if DIM == 2
            a_dt = min(dt[0], dt[1]);
        #endif
        #if DIM == 3 
            double a_dt_temp = min(dt[0], dt[1]);
            a_dt = min(a_dt_temp, dt[2]);
        #endif
	}

}
