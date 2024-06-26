#include <cstdio>
#include <cstring>
#include <cassert>
#include <cmath>

#include <vector>
#include <memory>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "Proto.H"
#include "SRK4.H"
#include "deformation.H"
#include <deform.hpp>
#include "Proto_WriteBoxData.H"
#include "Proto_Timer.H"
#include <array>
#include "writers.h"
#include <complex>
#include "Proto_Timer.H"
#define NUMCOMPS DIM+2

using namespace std;
using namespace Proto;



/***/
int main(int argc, char* argv[])
{
  PR_TIMER_SETFILE("proto.time.table");
  PR_TIMERS("main");
  int nx =128,ny = 128;
  int nx2 = 256,ny2 = 256;
  double dt = 0.0421875, dt2 = dt/2.0;
  double time = 0, timef = 0, time2 = 0;
  double z,z1,z2, dudx, dvdx, dvdy, dudy;
  double c = 0.75;
  double sumU = 0, sumVort = 0, sumVortm = 0, sumUg = 0;
  double tstop = 15; //15
  double maxStep,L1_eigen, L2_eigen;
  double R = 0.1125;
  double offset = 0.15;
  double h,hp,h_second,hp_second,x, y, x1,x2, ux, uy, ux0, uy0;
  double L = 1.5;
  double Amp = 4;
  vector<particle> X, X2;
  vector<double> t, angleError,lambda,lambda2, t_remap;
  vector<int> num_p;
  vector<double>  initial_radi;
  vector<double> v_theta,alpha, v_theta0;
  vector<array<double,DIM+12> > particles;
  vector<array<double, DIM+2>> corr;
  vector<array<double, DIM+9>> angle;
  vector<array<double, DIM>> corr2;
  vector<particle> Yp;
  vector<double> u,v, angleQ, eigenR, u2,v2, angleQ2, eigenR2,detF;
  vector<deform> Q, R2, vecFerror,Q2, R22, vecFerror2;
  vector<double> errorU,errorU1,errorVg, errorVg1, errorVort1,errorVort,errorVortm, errorpsi, errorLambda;
  vector<double> errorVortdt, errorVortpgp, errorDiffVg, errorDiffUg, errorDiffVort, errorDiffVortdt;  
  vector<double> errorpsi2,errorVg2, errorDiffV, errorDiffVdt, errorDiffV1, errorDiffV2,errorL2eigen, errorL1eigen;
  Point r;
  particle p, py;

  string solnfile, solnfilef, soln_prefix = "Vort2D";
  maxStep = tstop/dt;
  
  Box Bg(Point::Zeros(), Point({nx-1, ny-1}));
  Box Bg2(Point::Zeros(), Point({nx2-1, ny2-1}));
  BoxData<double> Ug(Bg);
  BoxData<double> Vg(Bg);
  BoxData<double> Ui(Bg);
  BoxData<double> Vi(Bg);
  BoxData<double> vort(Bg); BoxData<double> Lambda(Bg); BoxData<double> Lambda2(Bg2);
  BoxData<double> Vortg(Bg); BoxData<double> Vortg2(Bg2);
  BoxData<double> psi(Bg);
  Box Bp(Point::Zeros(), Point({4*(nx-1), 4*(nx-1)}));
  Box Bp_second(Point::Zeros(), Point({4*(nx2-1), 4*(nx2-1)}));
  string q,e;
  double Np = 0, wp = 0, Np2 = 0, wp2 = 0;;
  int M = log2(nx);
  interp interpol;
 
  string testname="test";
  string filen = "corr";
  string testname3 = "angle";
  string testname4 = "angle_omega_corr";

  const char *const varnames[] = {"x_1", "x_2", 
                                  "vorticity", 
				  "absvorterror", "relvorterror", 
				  "vortpgp", "vortpgp_error", "vortpgpdt", "vortpgpdt_error", "angleQ", "eigenR",
				   "detF", "lambdapgp", "lambdapgp_error"}; 

  h = L/(nx-1);
  hp = L/( 4*(nx-1) );
  h_second = L/(nx2 - 1);
  hp_second = L/( 4*(nx2-1) );

  double h2 = L/nx;
  
  Ug.setToZero(); Vg.setToZero();
  Vortg.setToZero(); vort.setToZero();
  psi.setToZero(); Vortg2.setToZero();
  Lambda.setToZero();

  //Initialize Position and Vorticity
  for( auto it = Bp.begin(); it != Bp.end(); it++){

    PR_TIMERS("main::deposition");

    r = it.operator*();
    p.x[0] = r[0]*hp-c;
    p.x[1] = r[1]*hp-c;

    z1 =( pow( pow((p.x[0]-offset ),2.0) + pow((p.x[1] ),2.0), 0.5 ))/R;
    z2 =( pow( pow((p.x[0]+offset ),2.0) + pow((p.x[1] ),2.0), 0.5 ))/R;

    if( z1 <= 1){
      p.strength = Amp*pow( ( 1 - pow(z1, 2.0) ), 7.0 );

      if( abs(p.strength) > pow(10, -6.0)){
        Np++;
        wp += p.strength*pow(hp, 2.0);
        lambda.push_back(1);
        X.push_back(p);

      }
    
    }

    if( z2 <= 1){
      p.strength = Amp*pow( ( 1 - pow(z2, 2.0) ), 7.0 );

      if( abs(p.strength) > pow(10, -7.0)){
        Np++;
        wp += p.strength*pow(hp, 2.0);
        lambda.push_back(1);
        X.push_back(p);

      }

    }



  }

    //Initialize Position and Vorticity
  for( auto it = Bp_second.begin(); it != Bp_second.end(); it++){


    r = it.operator*();
    p.x[0] = r[0]*hp_second-c;
    p.x[1] = r[1]*hp_second-c;

    z1 =( pow( pow((p.x[0]-offset ),2.0) + pow((p.x[1] ),2.0), 0.5 ))/R;
    z2 =( pow( pow((p.x[0]+offset ),2.0) + pow((p.x[1] ),2.0), 0.5 ))/R;

    if( z1 <= 1){
      p.strength = Amp*pow( ( 1 - pow(z1, 2.0) ), 7.0 );

      if( abs(p.strength) > pow(10, -6.0)){
        Np2++;
        wp2 += p.strength*pow(hp_second, 2.0);
        lambda2.push_back(1);
        X2.push_back(p);

      }

    }

    if( z2 <= 1){
      p.strength = Amp*pow( ( 1 - pow(z2, 2.0) ), 7.0 );

      if( abs(p.strength) > pow(10, -6.0)){
        Np2++;
        wp2 += p.strength*pow(hp_second, 2.0);
        lambda2.push_back(1);
        X2.push_back(p);

      }

    }



  }

  cout << " Np " << Np << " Np2 " << Np2 << endl;
  cout << " hp " << hp << " hp2 " << hp_second << endl;
  cout << " nx2 " << nx2 << endl;
  //Keep track of # particles, used to see growth in particles due to remapping  
  num_p.push_back(Np);

  State state(X,dt,hp, h, L, wp, lambda);
  deformation f(X); 
  RK4<State,F,DX> rk4;
  RK4<deformation, RHS, DF> rk4_f;
  State state_second(X2,dt2,hp_second, h_second, L, wp2, lambda2);
  deformation f_second(X2);

  int step = 0.16875/dt;
  int tpt;

  double sumU1 = 0, sumVort1 = 0, sumUg1 = 0;
  sumU = 0; sumVort = 0; sumUg = 0;
  Vortg.setToZero();

  {
  state.getVelocity(Lambda,Vortg,errorVg, errorpsi,u,v,Np,0,state);
  state_second.getVelocity(Lambda2,Vortg2,errorVg2, errorpsi2,u2,v2,Np2,0,state_second);

  }



  t_remap.push_back(0.0);
  array<double, DIM+12> temp;


  //Main time loop
  for( int k = 0;(k < maxStep);k++){
      PR_TIMERS("main::time_loop");


	 {

         //QR of deformation
         f.QR(Q,R2, angleQ, eigenR,detF);
        // f_second.QR(Q2,R22, angleQ2, eigenR2);
	 }
         vecFerror.clear(); angleError.clear(); particles.clear(); corr.clear();
         angle.clear(); corr2.clear();
         L1_eigen = 0; L2_eigen = 0;
        
	 //Calculate quantities to be written to particle reader
	 for( int i = 0; i <  state.X.size(); i++){

	    temp[0] = state.X.at(i).x[0]; temp[1] = state.X.at(i).x[1];
/*	    temp[2] = angleQ.at(i); temp[3] =eigenR.at(i);
	    temp[4] = f.F.at(i).F[0][0]; temp[5] = f.F.at(i).F[0][1];
            temp[6] = f.F.at(i).F[1][0]; temp[7] = f.F.at(i).F[1][1]; */
	    temp[2] = state.X.at(i).strength; temp[3] = state.corrVort_error.at(i);
	   
	    temp[4] = state.corrRelVort.at(i); 
	    temp[5] = state.vortpgp.at(i);
	    temp[6] = state.vortpgp_error.at(i);
	    temp[7] = state.vortpgpdt.at(i);
	    temp[8] = state.vortpgpdt_error.at(i);
            temp[9] = angleQ.at(i);
	    temp[10] = eigenR.at(i);
	    temp[11] = detF.at(i);
	    temp[12] = state.lambdapgp.at(i);
            temp[13] = state.lambdapgp_error.at(i);
	    particles.push_back(temp);
            L1_eigen += state.X.at(i).strength/Amp*abs(eigenR.at(i) - 1)*hp;
	    L2_eigen += pow( state.X.at(i).strength/Amp*abs(eigenR.at(i) - 1)*hp, 2.0  );
 

           }

	   L2_eigen = pow(L2_eigen, 0.5);
           
         {

         PWrite<DIM+12>(particles,varnames,testname,k);
         }  

        if(  L1_eigen > 0.1){     //(k+1)%10 == 0){
          state.remap();
//	  state_second.remap();
          num_p.push_back(state.X.size() );
          t_remap.push_back(time);
	  f.remap(state);
        }

//CURRENTLY NOT CALCULATING RATE
//         state.rate(state_second, time/dt);
	 {

	  errorVortpgp.push_back(state.sumVortpgp);
	  errorVortdt.push_back(state.sumVortpgpdt);
          errorDiffVort.push_back(state.diffVort);
          errorDiffVg.push_back(state.diffVg);
          errorDiffUg.push_back(state.diffUg);
	  errorDiffV.push_back(state.diffV);
          errorDiffVdt.push_back(state.diffVdt);
          errorDiffVortdt.push_back(state.diffVortdt);
	  errorDiffV1.push_back(state.diffV1);
	  errorDiffV2.push_back(state.diffV2);
	  errorLambda.push_back(state.sumLambdapgp);
	  errorL2eigen.push_back(L2_eigen);
	  errorL1eigen.push_back(L1_eigen);
         }

          t.push_back(time);
   
	  timef += dt;
          time += dt;



          {
          PR_TIMERS("main::RK4");

	  //NOT CALCULATING 2ND RATE
          //for(int i=0; i < 2;i++){
	  //   time2 += dt/2;
          //   rk4.advance(time2,dt2,state_second);
	  //   rk4_f.advance(time2,dt2,f_second); 
	  //}

	  
          rk4.advance(time,dt,state);
	  }
          

	  {
          PR_TIMERS("main::update_f");
	  //f_second.update_X(state_second);
	  f.update_X(state);

	  }



          rk4_f.advance(time,dt,f);

          
    }

  

  //Write curve files
  string errorfile2 = "errorVorticity.curve";
  ofstream f4(errorfile2);
  string errorfile4 = "errorVorticitym.curve";
  string errorfile5 = "errorVorticity1.curve";
  string errorfile7 = "particles.curve";
  ofstream f6(errorfile4);
  ofstream f7(errorfile5);
  ofstream f9(errorfile7);
  string eVortpgp = "errorVortpgp.curve";
  string eVortdt  = "errorVortdt.curve";
  ofstream f10(eVortpgp);
  ofstream f11(eVortdt);
  string  eDiffVort = "errorDiffVort.curve";
  string  eDiffVg   = "errorDiffVg.curve";
  string  eDiffUg   = "errorDiffUg.curve";
  string  eDiffV    = "errorDiffV.curve";
  string  eDiffVdt  = "errorDiffVdt.curve";
  string  eDiffVortdt  = "errorDiffVortdt.curve";
  string  eDiffV1      = "errorDiffV1.curve";
  string  eDiffV2      = "errorDiffV2.curve";
  string  lambdaerror      = "lambdaerror.curve";
  string  L2EigenError = "L2EigenError";
  string  L1EigenError = "L1EigenError";
  ofstream f12(eDiffVort);
  ofstream f13(eDiffVg);
  ofstream f14(eDiffUg);
  ofstream f15(eDiffV);
  ofstream f16(eDiffVdt);
  ofstream f17(eDiffVortdt);
  ofstream f18(eDiffV1);
  ofstream f19(eDiffV2);
  ofstream f20(lambdaerror);
  ofstream f21(L2EigenError);
  ofstream f22(L1EigenError);

  q = "# TIME";
  e = "# Vorticitypgp";


  if( f10.is_open() ) {

      //format in scientific notation
       f10.setf(std::ios::scientific, std::ios::floatfield);
       f10.precision(4);

       f10 << q << endl;


       f10 << e << endl;


      for(unsigned int i = 0; i < t.size(); i++){

          f10 << t.at(i) << " " << errorVortpgp.at(i)  <<  endl;
      }

   }

  e = "# Vorticitydt";


  if( f11.is_open() ) {

      //format in scientific notation
       f11.setf(std::ios::scientific, std::ios::floatfield);
       f11.precision(4);


      f11 << q << endl;


      f11 << e << endl;


      for(unsigned int i = 0; i < t.size(); i++){

          f11 << t.at(i) << " " << errorVortdt.at(i)  <<  endl;
      }

   }


  e = "# Particles";
   if( f9.is_open() ) {

      //format in scientific notation
       f9.setf(std::ios::scientific, std::ios::floatfield);
       f9.precision(4);


      f9 << q << endl;


      f9 << e << endl;

      for(unsigned int i = 0; i < t_remap.size(); i++){

          f9 << t_remap.at(i) << " " << num_p.at(i)  <<  endl;
      }

   }

  

   e = "# diffVort";


  if( f12.is_open() ) {

      //format in scientific notation
       f12.setf(std::ios::scientific, std::ios::floatfield);
       f12.precision(4);


      f12 << q << endl;


      f12 << e << endl;


      for(unsigned int i = 0; i < t.size(); i++){

          f12 << t.at(i) << " " << errorDiffVort.at(i)  <<  endl;
      }

   }

  cout << " diffvort" << endl;
   e = "# DiffVg";


   if( f13.is_open() ) {

      //format in scientific notation
       f13.setf(std::ios::scientific, std::ios::floatfield);
       f13.precision(4);


      f13 << q << endl;


      f13 << e << endl;


      for(unsigned int i = 0; i < t.size(); i++){

          f13 << t.at(i) << " " << errorDiffVg.at(i)  <<  endl;
      }

   }

   cout << "diffVg" << endl;
    e = "# DiffUg";


   if( f14.is_open() ) {

      //format in scientific notation
       f13.setf(std::ios::scientific, std::ios::floatfield);
       f13.precision(4);


      f14 << q << endl;


      f14 << e << endl;


      for(unsigned int i = 0; i < t.size(); i++){

          f14 << t.at(i) << " " << errorDiffUg.at(i)  <<  endl;
      }

   }

   e = "# DiffV";

   if( f15.is_open() ) {

      //format in scientific notation
       f15.setf(std::ios::scientific, std::ios::floatfield);
       f15.precision(4);


      f15 << q << endl;


      f15 << e << endl;


      for(unsigned int i = 0; i < t.size(); i++){

          f15 << t.at(i) << " " << errorDiffV.at(i)  <<  endl;
      }

   }


   e = "# DiffVdt";


   if( f16.is_open() ) {

      //format in scientific notation
       f16.setf(std::ios::scientific, std::ios::floatfield);
       f16.precision(4);


      f16 << q << endl;


      f16 << e << endl;


      for(unsigned int i = 0; i < t.size(); i++){

          f16 << t.at(i) << " " << errorDiffVdt.at(i)  <<  endl;
      }

   }

   
      e = "# DiffVortdt";


   if( f17.is_open() ) {

      //format in scientific notation
       f17.setf(std::ios::scientific, std::ios::floatfield);
       f17.precision(4);


      f17 << q << endl;


      f17 << e << endl;


      for(unsigned int i = 0; i < t.size(); i++){

          f17 << t.at(i) << " " << errorDiffVortdt.at(i)  <<  endl;
      }

   }

    e = "# DiffV1";


   if( f18.is_open() ) {

      //format in scientific notation
       f18.setf(std::ios::scientific, std::ios::floatfield);
       f18.precision(4);


      f18 << q << endl;


      f18 << e << endl;


      for(unsigned int i = 0; i < t.size(); i++){

          f18 << t.at(i) << " " << errorDiffV1.at(i)  <<  endl;
      }

   }

   e = "# DiffV2";

   if( f19.is_open() ) {

      //format in scientific notation
       f19.setf(std::ios::scientific, std::ios::floatfield);
       f19.precision(4);


      f19 << q << endl;


      f19 << e << endl;


      for(unsigned int i = 0; i < t.size(); i++){

          f19 << t.at(i) << " " << errorDiffV2.at(i)  <<  endl;
      }

   }

   e = "# LambdaError";

   if( f20.is_open() ) {

      //format in scientific notation
       f20.setf(std::ios::scientific, std::ios::floatfield);
       f20.precision(4);


      f20 << q << endl;


      f20 << e << endl;


      for(unsigned int i = 0; i < t.size(); i++){

          f20 << t.at(i) << " " << errorLambda.at(i)  <<  endl;
      }

   }


     e = "# L2EigenError";

   if( f21.is_open() ) {

      //format in scientific notation
       f21.setf(std::ios::scientific, std::ios::floatfield);
       f21.precision(4);


      f21 << q << endl;


      f21 << e << endl;


      for(unsigned int i = 0; i < t.size(); i++){

          f21 << t.at(i) << " " << errorL2eigen.at(i)  <<  endl;
      }

   }

    e = "# L1EigenError";

   if( f22.is_open() ) {

      //format in scientific notation
       f22.setf(std::ios::scientific, std::ios::floatfield);
       f22.precision(4);


      f22 << q << endl;


      f22 << e << endl;


      for(unsigned int i = 0; i < t.size(); i++){

          f22 << t.at(i) << " " << errorL1eigen.at(i)  <<  endl;
      }

   }



  cout << "diffUg" << endl;

PR_TIMER_REPORT();
}
