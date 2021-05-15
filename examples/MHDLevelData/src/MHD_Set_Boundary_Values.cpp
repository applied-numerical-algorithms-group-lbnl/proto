#include "Proto.H"
#include "MHD_Set_Boundary_Values.H"
#include "CommonTemplates.H"
#include "Proto_Timer.H"
#include "Proto_WriteBoxData.H"
#include "Proto_LevelBoxData.H"
#include "Proto_ProblemDomain.H"
#include "MHDOp.H"
#include "MHD_Mapping.H"

typedef BoxData<double,1,1,1> Scalar;
typedef BoxData<double,NUMCOMPS,1,1> Vector;

namespace MHD_Set_Boundary_Values {

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
	void fill_ghostsF(State& a_U,
	                  V& a_x,
	                  const double a_gamma)
	{
		double x = a_x(0);
		double y = a_x(1);
		double gamma = a_gamma;
		double rho = 0.0;
		double p = 0.0;
		double u = 0.0;
		double v = 0.0;
		double w = 0.0;
		double Bx = 0.0;
		double By = 0.0;
		double Bz = 0.0;



		double theta = atan2(y,x);
		// rho = 1.0;
		// p = 1.0;
		// u = 1.0;
		// v = 1.0;
		double rho_0 = 1.4;
		double delta_rho_0 = 0.14;
		rho = rho_0 + delta_rho_0*exp(-400.0*(x-0.5)*(x-0.5))*pow(cos(PI*(x-0.5)),16.0);
		p = pow(rho/rho_0,gamma);

		double e = p/(gamma-1.0) + rho*(u*u+v*v+w*w)/2.0 + (Bx*Bx+By*By+Bz*Bz)/8.0/PI;
		//cout << "here" << endl;
#if DIM == 2
		a_U(0) = rho;      //rho
		a_U(1) = rho*u;    //Momentum-x
		a_U(2) = rho*v;    //Momentum-y
		a_U(3) = e;        //Energy
		a_U(4) = Bx;       //Bx
		a_U(5) = By;       //By
#endif

#if DIM == 3
		a_U(0) = rho;      //rho
		a_U(1) = rho*u;    //Momentum-x
		a_U(2) = rho*v;    //Momentum-y
		a_U(3) = rho*w;    //Momentum-z
		a_U(4) = e;        //Energy
		a_U(5) = Bx;       //Bx
		a_U(6) = By;       //By
		a_U(7) = Bz;       //Bz
#endif
	}
	PROTO_KERNEL_END(fill_ghostsF, fill_ghosts)


	PROTO_KERNEL_START
	void iotaFuncF(Point           & a_p,
	               V               & a_X,
	               double a_h)
	{
		for (int ii = 0; ii < DIM; ii++)
		{
			a_X(ii) = a_p[ii]*a_h + 0.5*a_h;
		}
	}
	PROTO_KERNEL_END(iotaFuncF,iotaFunc)




	void Set_Boundary_Values(BoxData<double,NUMCOMPS>& a_JU,
	                         const Box& a_rangeBox,
	                         const ProblemDomain& a_probDomain,
	                         const double a_dx,
	                         const double a_gamma,
	                         const BoxData<double,1>& Jacobian_ave,
	                         const int a_lowBCtype,
	                         const int a_highBCtype)
	{

#if DIM == 2
		Box dbx0 = a_JU.box();
		//Filling ghost cells for low side of dir == 0 here. This will be the inner boundary in r direction once we map to polar, spherical, or cubed sphere grids.
		if (dbx0.low()[0] < a_probDomain.box().low()[0] && a_lowBCtype == 1) {
			Point ghost_low = dbx0.low();
			Point ghost_high = Point(a_probDomain.box().low()[0]-1, dbx0.high()[1]);
			Box BoundBox(ghost_low,ghost_high);

			//Box BoundBox1 = BoundBox.grow(1);
			Vector a_U(BoundBox);
			BoxData<double,DIM> eta(BoundBox);
			forallInPlace_p(iotaFunc, BoundBox, eta, a_dx);
			BoxData<double,DIM> x(BoundBox);
			MHD_Mapping::eta_to_x_calc(x,eta);
			forallInPlace(fill_ghosts,BoundBox,a_U,x,a_gamma);
			//Scalar Jacobian_ave(BoundBox1);
			//MHD_Mapping::Jacobian_Ave_calc(Jacobian_ave,a_dx,BoundBox1);
			Vector a_JU_ghost = forall<double,NUMCOMPS>(dot_pro_calcF, Jacobian_ave, a_U);
			a_JU_ghost.copyTo(a_JU,BoundBox);
		}


		//Filling ghost cells for high side of dir == 0 here. This will be the outer boundary in r direction once we map to polar, spherical, or cubed sphere grids.
		if (dbx0.high()[0] > a_probDomain.box().high()[0] && a_highBCtype == 1) {
			Point ghost_low = Point(a_probDomain.box().high()[0]+1, dbx0.low()[1]);
			Point ghost_high = dbx0.high();
			Box BoundBox(ghost_low,ghost_high);
			//cout << BoundBox.low()[0] << " " <<  BoundBox.low()[1] << " " <<  BoundBox.high()[0] << " " <<  BoundBox.high()[1] << endl;
			//Box BoundBox1 = BoundBox.grow(1);
			Vector a_U(BoundBox);
			BoxData<double,DIM> eta(BoundBox);
			forallInPlace_p(iotaFunc, BoundBox, eta, a_dx);
			BoxData<double,DIM> x(BoundBox);
			MHD_Mapping::eta_to_x_calc(x,eta);
			forallInPlace(fill_ghosts,BoundBox,a_U,x,a_gamma);
			//Scalar Jacobian_ave(BoundBox1);
			//MHD_Mapping::Jacobian_Ave_calc(Jacobian_ave,a_dx,BoundBox1);
			Vector a_JU_ghost = forall<double,NUMCOMPS>(dot_pro_calcF, Jacobian_ave, a_U);
			a_JU_ghost.copyTo(a_JU,BoundBox);
		}



		//Filling ghost cells for low side of dir == 0 here, to make an open boundary. Ghost cells are filled by innermost cells of problem domain. This will be the inner boundary in r direction once we map to polar, spherical, or cubed sphere grids.
		if (dbx0.low()[0] < a_probDomain.box().low()[0] && a_lowBCtype == 2) {
			Point source_low = Point(a_probDomain.box().low()[0],a_JU.box().low()[1]);
			Point source_high = Point(a_probDomain.box().low()[0],a_JU.box().high()[1]);
			Box sourceBox(source_low,source_high);
			Box sourceBox1 = sourceBox.grow(1);

			//Scalar Jacobian_ave(sourceBox1);
			//MHD_Mapping::Jacobian_Ave_calc(Jacobian_ave,a_dx,sourceBox1);
			Vector a_U(sourceBox1);
			MHD_Mapping::JU_to_U_2ndOrdercalc(a_U, a_JU, Jacobian_ave, sourceBox1);
			Point ghost_low = dbx0.low();
			Point ghost_high = Point(a_probDomain.box().low()[0]-1, dbx0.high()[1]);
			Box BoundBox(ghost_low,ghost_high);
			Box BoundBox1 = BoundBox.grow(1);
			Vector a_U_ghost(BoundBox);
			Scalar Jacobian_ave2(BoundBox1);
			MHD_Mapping::Jacobian_Ave_calc(Jacobian_ave2,a_dx,BoundBox1);

			for (int i = 1; i <= NGHOST; i++ ) {
				a_U.copyTo(a_U_ghost,sourceBox,Point::Basis(0)*(-i));// Using shifting option of copyTo
			}
			Vector a_JU_ghost = forall<double,NUMCOMPS>(dot_pro_calcF, Jacobian_ave2, a_U_ghost);
			a_JU_ghost.copyTo(a_JU,BoundBox);
		}





		//Filling ghost cells for high side of dir == 0 here, to make an open boundary. Ghost cells are filled by outermost cells of problem domain. This will be the outer boundary in r direction once we map to polar, spherical, or cubed sphere grids.
		if (dbx0.high()[0] > a_probDomain.box().high()[0] && a_highBCtype == 2) {
			Point source_low = Point(a_probDomain.box().high()[0],a_JU.box().low()[1]);
			Point source_high = Point(a_probDomain.box().high()[0],a_JU.box().high()[1]);
			Box sourceBox(source_low,source_high);
			Box sourceBox1 = sourceBox.grow(1);

			//Scalar Jacobian_ave(sourceBox1);
			//MHD_Mapping::Jacobian_Ave_calc(Jacobian_ave,a_dx,sourceBox1);
			Vector a_U(sourceBox1);
			MHD_Mapping::JU_to_U_2ndOrdercalc(a_U, a_JU, Jacobian_ave, sourceBox1);
			Point ghost_low = Point(a_probDomain.box().high()[0]+1, dbx0.low()[1]);
			Point ghost_high = dbx0.high();
			Box BoundBox(ghost_low,ghost_high);
			Box BoundBox1 = BoundBox.grow(1);
			Vector a_U_ghost(BoundBox);
			Scalar Jacobian_ave2(BoundBox1);
			MHD_Mapping::Jacobian_Ave_calc(Jacobian_ave2,a_dx,BoundBox1);

			for (int i = 1; i <= NGHOST; i++ ) {
				a_U.copyTo(a_U_ghost,sourceBox,Point::Basis(0)*(i));// Using shifting option of copyTo
			}
			Vector a_JU_ghost = forall<double,NUMCOMPS>(dot_pro_calcF, Jacobian_ave2, a_U_ghost);
			a_JU_ghost.copyTo(a_JU,BoundBox);
		}

#endif
	}




	void Set_Jacobian_Values(BoxData<double,1>& a_Jacobian_ave,
	                         const Box& a_rangeBox,
	                         const ProblemDomain& a_probDomain,
	                         const double a_dx,
	                         const double a_gamma,
	                         const int a_lowBCtype,
	                         const int a_highBCtype)
	{

#if DIM == 2
		Box dbx0 = a_Jacobian_ave.box();
		//Filling Jacobian values in ghost cells for low side of dir == 0 here. This will be the inner boundary in r direction once we map to polar, spherical, or cubed sphere grids.
		if (dbx0.low()[0] < a_probDomain.box().low()[0] && a_lowBCtype != 0) {
			Point ghost_low = dbx0.low();
			Point ghost_high = Point(a_probDomain.box().low()[0]-1, dbx0.high()[1]);
			Box BoundBox(ghost_low,ghost_high);
			Box BoundBox1 = BoundBox.grow(1);
			Scalar Jacobian_ave2(BoundBox1);
			MHD_Mapping::Jacobian_Ave_calc(Jacobian_ave2,a_dx,BoundBox1);
			Jacobian_ave2.copyTo(a_Jacobian_ave,BoundBox);
		}


		//Filling Jacobian values in ghost cells for high side of dir == 0 here. This will be the outer boundary in r direction once we map to polar, spherical, or cubed sphere grids.
		if (dbx0.high()[0] > a_probDomain.box().high()[0] && a_highBCtype != 0) {
			Point ghost_low = Point(a_probDomain.box().high()[0]+1, dbx0.low()[1]);
			Point ghost_high = dbx0.high();
			Box BoundBox(ghost_low,ghost_high);
			Box BoundBox1 = BoundBox.grow(1);
			Scalar Jacobian_ave2(BoundBox1);
			MHD_Mapping::Jacobian_Ave_calc(Jacobian_ave2,a_dx,BoundBox1);
			Jacobian_ave2.copyTo(a_Jacobian_ave,BoundBox);
		}

		//std::string filename="Jacobian_ave";
		//WriteBoxData(filename.c_str(),a_Jacobian_ave,a_dx);
#endif
	}





}