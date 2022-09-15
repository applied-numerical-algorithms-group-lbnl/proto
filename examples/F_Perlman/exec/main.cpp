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
  int nx = 64,ny = 64;
  double dt = 0.084375;
  double time = 0, timef = 0;
  double z, z2, dudx, dvdx, dvdy, dudy;
  double c = 2.0;
  double sumU = 0, sumVort = 0, sumVortm = 0, sumUg = 0;
  double tstop = 20;
  double maxStep;
  double R = 1;
  double h,hp,x, y, x1,x2, ux, uy, ux0, uy0;
  double L = 4.0;

  vector<particle> X;
  vector<double> t, angleError,lambda, t_remap;
  vector<int> num_p;
  vector<double>  initial_radi;
  vector<double> v_theta,alpha, v_theta0;
  vector<array<double,DIM+17> > particles;
  vector<array<double, DIM+2>> corr;
  vector<array<double, DIM+9>> angle;
  vector<array<double, DIM>> corr2;
  vector<particle> Yp;
  vector<double> u,v, angleQ, eigenR;
  vector<deform> Q, R2, vecFerror;
  vector<double> errorU,errorU1,errorVg, errorVg1, errorVort1,errorVort,errorVortm, errorpsi;
  
  Point r;
  particle p, py;

  string solnfile, solnfilef, soln_prefix = "Vort2D";
  maxStep = tstop/dt;
  
  Box Bg(Point::Zeros(), Point({nx-1, ny-1}));
  BoxData<double> Ug(Bg);
  BoxData<double> Vg(Bg);
  BoxData<double> Ui(Bg);
  BoxData<double> Vi(Bg);
  BoxData<double> vort(Bg);
  BoxData<double> Vortg(Bg);
  BoxData<double> psi(Bg);
  Box Bp(Point::Zeros(), Point({4*(nx-1), 4*(nx-1)}));
  string q,e;
  double Np = 0, wp = 0;
  int M = log2(nx);
  interp interpol;
 
  const char *const varnames3[] = {"alpha1", "alpha2", "QTx1", "QTx2","error_1", "error_2", "alpha1_norm", "angle", "QTx_magn", "alpha_norm_magn", "error_magn"}; 
  string testname="test";
  string filen = "corr";
  string testname3 = "angle";
  string testname4 = "angle_omega_corr";

  const char *const varnames4[] = {"LambdaVort_error", "Vort_error"};
  const char *const varnames[] = {"x_1", "x_2", "angleQ", "eigenR", "F11",
                                  "F12", "F21", "F22", "vorticity", 
				  "absvorterror", "relvorterror", "lambda", "lambda_abserror", 
				  "QTx1", "position_error_x1", "position_error_x2", "position_error_norm", 
				  "alpha_norm1", "radi_error"}; 
  const char *const varnames2[] = {"angleQ", "eigenR", "absvorterror", "relvorterror"};

  h = L/(nx-1);
  hp = L/( 4*(nx-1) );
  double h2 = L/nx;
  
  Ug.setToZero(); Vg.setToZero();
  Vortg.setToZero(); vort.setToZero();
  psi.setToZero();


 //Calculate exact velocity on grid
 for( auto it = Bg.begin(); it != Bg.end(); it++){

    r = *it;

    p.x[0] = r[0]*h-c;
    p.x[1] = r[1]*h-c;
    z = ( pow( pow((p.x[0] ),2.0) + pow((p.x[1] ),2.0), 0.5 ));

    if( z <= R){
        vort(r) +=  pow( ( 1 - pow(z, 2.0) ), 7.0 );
        Ug(r) += -1.0/(16.0*pow(z, 2.0) )*(1 - pow( (1 - pow(z, 2.0) ), 8.0) )*(p.x[1] );
        Vg(r) += 1.0/(16.0*pow(z, 2.0) )*(1 - pow( (1 - pow(z, 2.0) ), 8.0) )*(p.x[0] );

    } else{

         z2 = pow((p.x[0] ),2.0) + pow((p.x[1] ),2.0);
         Ug(r) += -1.0/(16.0*z2) *(p.x[1] );
         Vg(r) +=  1.0/(16.0*z2) *(p.x[0] );
    }

  }


  //Initialize Position and Vorticity
  for( auto it = Bp.begin(); it != Bp.end(); it++){

    PR_TIMERS("main::deposition");

    r = it.operator*();
    p.x[0] = r[0]*hp-c;
    p.x[1] = r[1]*hp-c;
    z = ( pow( pow((p.x[0] ),2.0) + pow((p.x[1] ),2.0), 0.5 ));


    if( z <= 1){
      p.strength = pow( ( 1 - pow(z, 2.0) ), 7.0 );

      if( p.strength > pow(10, -7.0)){
        Np++;
        initial_radi.push_back(z);
        wp += p.strength*pow(hp, 2.0);
        lambda.push_back(1);
        X.push_back(p);

      }
    
    }


  }

  //Keep track of # particles, used to see growth in particles due to remapping  
  num_p.push_back(Np);

  #ifdef PR_HDF5
     HDF5Handler h5;
     h5.writePatch(h,Ug,{"U"}, "exact_U%i.hdf5", 0);
     h5.writePatch(h,Vg,{"V"}, "exact_V%i.hdf5", 0);
  #endif

  State state(X,dt,hp, h, L, wp, lambda);
  deformation f(X); 
  RK4<State,F,DX> rk4;
  RK4<deformation, RHS, DF> rk4_f;
  int step = 0.16875/dt;
  int tpt;

  double sumU1 = 0, sumVort1 = 0, sumUg1 = 0;
  sumU = 0; sumVort = 0; sumUg = 0;
  Vortg.setToZero();

  {
  state.getVelocity(Vortg,errorVg, errorpsi,u,v,Np,0,state);
  }

  sumVortm = 0;
  sumVort = 0; sumVort1 = 0;
  for(auto it = Bg.begin(); it != Bg.end(); it++){
     
      r = it.operator*();
      sumVort +=  pow( *Vortg.data(r) - *vort.data(r), 2.0)*pow(h, 2.0);
      sumVort1 +=  abs( *Vortg.data(r) - *vort.data(r))*pow(h, 2.0);
      if( abs( *Vortg.data(r) - *vort.data(r)) > sumVortm){
        sumVortm = abs( *Vortg.data(r) - *vort.data(r));

      }

  }

  sumU = pow(sumU, 0.5);
  sumVort = pow(sumVort, 0.5);
  errorU.push_back(sumU); errorVort.push_back(sumVort);
  errorU1.push_back(sumU1); errorVort1.push_back(sumVort1);
  errorVg.push_back(state.sumUg); errorVg1.push_back(state.sumUg1);
  errorVortm.push_back(sumVortm);

  t.push_back(0); t_remap.push_back(0.0);
  array<double, DIM+17> temp;
  array<double, DIM+2> temp2;
  array<double, DIM+9> temp3;
  array<double, DIM> temp4;


  //Main time loop
  for( int k = 0;(k < maxStep);k++){
      PR_TIMERS("main::time_loop");


	 {
         //QR of deformation
         f.QR(Q,R2, angleQ, eigenR);
	 }
         vecFerror.clear(); angleError.clear(); particles.clear(); corr.clear();
         angle.clear(); corr2.clear();
    

	 //Calculate quantities to be written to particle reader
	 for( int i = 0; i <  state.X.size(); i++){


            z2 = pow( state.X.at(i).x[0], 2.0) + pow(state.X.at(i).x[1], 2.0);   

	    x = state.X.at(i).x[0]; y = state.X.at(i).x[1];
            z = ( pow( pow((x ),2.0) + pow((y ),2.0), 0.5 ));

	    if( z <= R){
            v_theta.push_back( -1.0/(16.0*z2)*(1 - pow( (1 - z2), 8.0) ));

            } else{
 

            v_theta.push_back(( -1.0/(16.0*z2) ));

            }
             
	    alpha.push_back( pow( pow( X.at(i).x[0], 2.0) + pow( X.at(i).x[1], 2.0), 0.5));

 
	    if( z2 > pow( 10.0, -10.0) ){

	    x1 = cos(time*v_theta.at(i))*(state.X.at(i).x[0]) + sin(time*v_theta.at(i))*(state.X.at(i).x[1]);
	    x2 = cos(time*v_theta.at(i))*(state.X.at(i).x[1]) - sin(time*v_theta.at(i))*(state.X.at(i).x[0]);

	    x1 = x1;
	 
	    x2 = x2;


	    temp[0] = state.X.at(i).x[0]; temp[1] = state.X.at(i).x[1];
	    temp[2] = angleQ.at(i); temp[3] =eigenR.at(i);
	    temp[4] = f.F.at(i).F[0][0]; temp[5] = f.F.at(i).F[0][1];
            temp[6] = f.F.at(i).F[1][0]; temp[7] = f.F.at(i).F[1][1];
	    temp[8] = state.X.at(i).strength; temp[9] = state.corrVort_error.at(i);
	   
	    temp[10] = state.corrRelVort.at(i); temp[11] = state.lambda.at(i); 
	    temp[12] = abs(abs(state.lambda.at(i) - 1)*state.X.at(i).strength);
	    temp[14] = (abs(x1) - abs(X.at(i).x[0]) );
	    temp[15] = (abs(x2) - abs(X.at(i).x[1]) );

            temp[16] = pow( pow( abs(x1) - abs(X.at(i).x[0]), 2.0) + pow( abs(x2) - abs(X.at(i).x[1]), 2.0), 0.5);
	    temp[17] = state.X.at(i).x[0]/alpha.at(i);
	    temp[18] = abs( z-initial_radi.at(i));
	    temp[13] = x1;
            temp2[0] = angleQ.at(i); temp2[1] = eigenR.at(i);
	    temp2[2] = state.corrVort_error.at(i);
	    temp2[3] = state.corrRelVort.at(i);

	    temp3[0] = X.at(i).x[0];
	    temp3[1] = X.at(i).x[1];
	    temp3[2] = x1;
	    temp3[3] = x2;
	    temp3[4] = (abs(x1) - abs(X.at(i).x[0]) );
	    temp3[5] = (abs(x2) - abs(X.at(i).x[1]) ); 
	    temp3[6] = X.at(i).x[0]/alpha.at(i);
	    temp3[7] = abs(v_theta.at(i)*time); // X.at(i).x[1]/alpha.at(i);
            temp3[8] = pow( pow(x1, 2.0) + pow(x2, 2.0), 0.5);
	    temp3[9] = pow( pow(X.at(i).x[0]/alpha.at(i), 2.0) + pow(X.at(i).x[1]/alpha.at(i), 2.0), 0.5);
	    temp3[10] = pow( pow( abs(x1) - abs(X.at(i).x[0]/alpha.at(i)), 2.0) + pow( abs(x2) - abs(X.at(i).x[1]/alpha.at(i)), 2.0), 0.5);
	  
	    temp4[0] = abs(abs(state.lambda.at(i) - 1)*state.X.at(i).strength);
	    temp4[1] = state.corrVort_error.at(i);
	  
	    if( temp2[2] > pow(10, -3.0)){
	      corr.push_back(temp2);
	    }

	    corr2.push_back(temp4);

	    particles.push_back(temp);
	    angle.push_back(temp3);
 
	    corr2.push_back(temp4);

           }

	 }	 

         {

         PWrite<DIM+17>(particles,varnames,testname,k);
         PWrite<DIM+2>(corr, varnames2, filen, k);
	 PWrite<DIM+9>(angle, varnames3, testname3, k);
         PWrite<DIM>(corr2, varnames4, testname4, k);
         }  

	 {
          errorU.push_back(state.sumU);
          errorVort.push_back(state.sumVort);
          errorVg.push_back(state.sumUg );
          errorU1.push_back(state.sumU1);
          errorVort1.push_back(state.sumVort1);
          errorVg1.push_back(state.sumUg1);
          errorVortm.push_back(state.sumVortm);
         }

          t.push_back(time);
   
	  timef += dt;
          time += dt;


	//if( (k+1)%6 == 0){
        //  state.remap();
        //  f.remap( state);
        //  timef = 0;
	//  num_p.push_back(state.X.size() );
	//  t_remap.push_back(time);
        //}

          {
          PR_TIMERS("main::RK4");

          rk4.advance(time,dt,state);
	  }


	  {
          PR_TIMERS("main::update_f");
	  f.update_X(state);
	  }

	  {
	   PR_TIMERS("main::RK4_f");
	   rk4_f.advance(time, dt, f);
	  }


	  

          
    }

  {

  //Write curve files
  string errorfile = "errorVelocity.curve";
  string errorfile2 = "errorVorticity.curve";
  string errorfile3 = "errorVg.curve";
  ofstream f3(errorfile);
  ofstream f4(errorfile2);
  ofstream f5(errorfile3);
  string errorfile4 = "errorVorticitym.curve";
  string errorfile5 = "errorVorticity1.curve";
  string errorfile6 = "errorVg1.curve";
  string errorfile7 = "particles.curve";
  ofstream f6(errorfile4);
  ofstream f7(errorfile5);
  ofstream f8(errorfile6);
  ofstream f9(errorfile7);

  q = "# TIME";
  e = "# Velocity";
 
  if( f3.is_open() ) {

      //format in scientific notation
       f3.setf(std::ios::scientific, std::ios::floatfield);
       f3.precision(4);

       f3 << q << endl;


       f3 << e << endl;


      for(unsigned int i = 0; i < t.size(); i++){

          f3 << t.at(i) << " " << errorU.at(i)  <<  endl;

      }

   }


  e = "# Vorticity";


  if( f4.is_open() ) {

      //format in scientific notation
       f4.setf(std::ios::scientific, std::ios::floatfield);
       f4.precision(4);


      f4 << q << endl;


      f4 << e << endl;


      for(unsigned int i = 0; i < t.size(); i++){

          f4 << t.at(i) << " " << errorVort.at(i)  <<  endl;

      }

   }

   e = "# Velocityg";


  if( f5.is_open() ) {

      //format in scientific notation
       f5.setf(std::ios::scientific, std::ios::floatfield);
       f5.precision(4);


      f5 << q << endl;


      f5 << e << endl;


      for(unsigned int i = 0; i < t.size(); i++){

          f5 << t.at(i) << " " << errorVg.at(i)  <<  endl;

      }

   }


  e = "# Vorticitym";


  if( f6.is_open() ) {

      //format in scientific notation
       f6.setf(std::ios::scientific, std::ios::floatfield);
       f6.precision(4);


      f6 << q << endl;


      f6 << e << endl;


      for(unsigned int i = 0; i < t.size(); i++){

          f6 << t.at(i) << " " << errorVortm.at(i)  <<  endl;

      }

   }

  e = "# Vort1";


  if( f7.is_open() ) {

      //format in scientific notation
       f7.setf(std::ios::scientific, std::ios::floatfield);
       f7.precision(4);


      f7 << q << endl;


      f7 << e << endl;


      for(unsigned int i = 0; i < t.size(); i++){

          f7 << t.at(i) << " " << errorVort1.at(i)  <<  endl;

      }

   }

  e = "# Vg1";


  if( f8.is_open() ) {

      //format in scientific notation
       f8.setf(std::ios::scientific, std::ios::floatfield);
       f8.precision(4);


      f8 << q << endl;


      f8 << e << endl;


      for(unsigned int i = 0; i < t.size(); i++){

          f8 << t.at(i) << " " << errorVg1.at(i)  <<  endl;

      }

   }

  e = "# Particles";
   if( f9.is_open() ) {

      //format in scientific notation
       f9.setf(std::ios::scientific, std::ios::floatfield);
       f9.precision(4);


      f9 << q << endl;


      f9 << e << endl;
   }

  }
PR_TIMER_REPORT();
}
