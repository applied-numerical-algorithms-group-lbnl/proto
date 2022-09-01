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
  int nx = 32;
  int ny = 32;
  Box Bg(Point::Zeros(), Point({nx-1, ny-1}));
  vector<particle> X;
  double dt = 0.16875;
  double time = 0, timef = 0;
  vector<double> t, angleError, t_remap;
  vector<int> num_p;
  vector<array<double,DIM+11> > particles;
  vector<array<double, DIM+2>> corr;
  deform d;
  double tstop = 20;
  vector<double> omeg;
  double maxStep;
  Point r;
  double h;
  double R = 1;
  double hp, x, y;
  complex<double> beta,sqrtv, sqrtu, c_dudy,c_dudx, c_dvdy, c_dvdx, alpha,F12, F11, F21, F22, lambda1, lambda2;
  double L = 3;
  BoxData<double> Ug(Bg);
  BoxData<double> Vg(Bg);
  BoxData<double> Ui(Bg);
  BoxData<double> Vi(Bg);
  particle p;
  string solnfile;
  string solnfilef;
  string soln_prefix = "Vort2D";
  string solnf_prefix = "F";
  maxStep = tstop/dt;
  double z, z2, dudx, dvdx, dvdy, dudy;
  double c = 1.5;
  double sumU = 0;
  double sumVort = 0, sumVortm =0;
  double sumUg = 0;
  vector<double> u,v, angleQ, eigenR;
  vector<deform> Q, R2, vecFerror;
  deform errorF;
  BoxData<double> vort(Bg);
  BoxData<double> Vortg(Bg);
  BoxData<double> psi(Bg);
  Box Bp(Point::Zeros(), Point({4*(nx-1), 4*(nx-1)}));
  BoxData<double, DIM> F_result;
  vector<double> errorU,errorU1,errorVg, errorVg1, errorVort1,errorVort,errorVortm, errorpsi;
  string q,e;
  double Np = 0;
  double wp = 0;
  int M = log2(nx);
  
  string testname="test";
  string filen = "corr";
  const char *const varnames[] = {"x_1", "x_2", "angleQ", "eigenR", "F11",
                                  "F12", "F21",
                                  "F22", "vorticity", "absvorterror", "relvorterror", "omega", "omega_abserror"}; 
  const char *const varnames2[] = {"angleQ", "eigenR", "absvorterror", "relvorterror"};

  h = L/(nx-1);
  hp = L/( 4*(nx-1) );
  double h2 = L/nx;
  Hockney hockney(h, M);
  Ug.setToZero(); Vg.setToZero();
  Vortg.setToZero(); vort.setToZero();
  psi.setToZero();

 for( auto it = Bg.begin(); it != Bg.end(); it++){

    r = it.operator*();

    p.x[0] = r[0]*h;
    p.x[1] = r[1]*h;

    z = ( pow( pow((p.x[0] - c),2.0) + pow((p.x[1] - c),2.0), 0.5 ));

    if( z <= R){
        vort(r) +=  pow( ( 1 - pow(z, 2.0) ), 7.0 );

        Ug(r) += -1.0/(16.0*pow(z, 2.0) )*(1 - pow( (1 - pow(z, 2.0) ), 8.0) )*(p.x[1] - c);
        Vg(r) += 1.0/(16.0*pow(z, 2.0) )*(1 - pow( (1 - pow(z, 2.0) ), 8.0) )*(p.x[0] - c);



    } else{

         z2 = pow((p.x[0] - c),2.0) + pow((p.x[1] - c),2.0);

         Ug(r) += -1.0/(16.0*z2) *(p.x[1] - c);
         Vg(r) +=  1.0/(16.0*z2) *(p.x[0] - c);




    }

  }

 vort.copyTo(psi);
 {
  hockney.convolve(psi);
  PR_TIMERS("main::hockney");
 }

   //Initialize Position and Vorticity
  for( auto it = Bp.begin(); it != Bp.end(); it++){

    r = it.operator*();

    p.x[0] = r[0]*hp;
    p.x[1] = r[1]*hp;

    z = ( pow( pow((p.x[0] - c),2.0) + pow((p.x[1] - c),2.0), 0.5 ));


    if( z <= 1){
    
    	    
       p.strength = pow( ( 1 - pow(z, 2.0) ), 7.0 );

       if( p.strength > pow(10, -7.0)){
       Np++;

       wp += p.strength*pow(hp, 2.0);
       X.push_back(p);
       omeg.push_back(1);
       }
    }

    PR_TIMERS("main::deposition");
  }

num_p.push_back(Np);

#ifdef PR_HDF5
     HDF5Handler h5;
     h5.writePatch(h,Ug,{"U"}, "exact_U%i.hdf5", 0);
     h5.writePatch(h,Vg,{"V"}, "exact_V%i.hdf5", 0);

#endif

  State state(X,dt,hp, h, L, wp, omeg);
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
  PR_TIMERS("main::get_velocity");
  }

  for(int i = 0; i < Np; i++){
    z =( pow( pow((state.X.at(i).x[0] - c),2.0) + pow((state.X.at(i).x[1] - c),2.0), 0.5 ));

    if( z <= 1){

       sumU += ( pow( u[i]-(state.X.at(i).x[1]-c)/(16*pow(z, 2.0) )*(1-pow(1-pow(z, 2.0), 8.0 )), 2.0 ) + pow( v[i]+(state.X.at(i).x[0]-c)/(16*pow(z, 2.0) )*(1-pow(1-pow(z, 2.0), 8.0 )), 2.0 ))*pow(hp, 2.0);

       sumU1 += ( abs( u[i]-(state.X.at(i).x[1]-c)/(16*pow(z, 2.0) )*(1-pow(1-pow(z, 2.0), 8.0 )) ) + abs( v[i]+(state.X.at(i).x[0]-c)/(16*pow(z, 2.0) )*(1-pow(1-pow(z, 2.0), 8.0 )) ))*pow(hp, 2.0);

    }else{

       z2 = pow((p.x[0] - c),2.0) + pow((p.x[1] - c),2.0);

       sumU += ( pow( u[i]-(state.X.at(i).x[1]-c)/(16*pow(z, 2.0) ), 2.0 ) + pow( v[i]+(state.X.at(i).x[0]-c)/(16*pow(z, 2.0) ), 2.0 ) )*pow(hp, 2.0);
       sumU1 += ( abs( u[i]-(state.X.at(i).x[1]-c)/(16*pow(z, 2.0))) + abs( v[i]+(state.X.at(i).x[0]-c)/(16*pow(z, 2.0) )) )*pow(hp, 2.0);

    }



   }

  sumVortm =0;
  for(auto it = Bg.begin(); it != Bg.end(); it++){
      r = it.operator*();
      sumVort +=  pow( *Vortg.data(r) - *vort.data(r), 2.0)*pow(h, 2.0);
      sumVort1 +=  abs( *Vortg.data(r) - *vort.data(r))*pow(h, 2.0);
      PR_TIMERS("main::vorticity_error");
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
  array<double, DIM+11> temp;
  array<double, DIM+2> temp2;


  for( int k = 0;(k < maxStep);k++)
  {
	 {
         f.QR(Q,R2, angleQ, eigenR);
         PR_TIMERS("main::QR");
	 }
         vecFerror.clear(); angleError.clear(); particles.clear(); corr.clear();

       for( int i = 0; i <  state.X.size(); i++){

            z2 = pow( state.X.at(i).x[0]-1.5, 2.0) + pow(state.X.at(i).x[1]-1.5, 2.0);   

	    x = state.X.at(i).x[0]; y = state.X.at(i).x[1];
            z = ( pow( pow((x - c),2.0) + pow((y - c),2.0), 0.5 ));
    
	    if( z2 > pow( 10.0, -10.0) ){
	    temp[0] = state.X.at(i).x[0]; temp[1] = state.X.at(i).x[1];
	    temp[2] = angleQ.at(i); temp[3] =eigenR.at(i);
	    temp[4] = f.F.at(i).F[0][0]; temp[5] = f.F.at(i).F[0][1];
            temp[6] = f.F.at(i).F[1][0]; temp[7] = f.F.at(i).F[1][1];
	    temp[8] = state.X.at(i).strength; temp[9] = state.corrVort_error.at(i);
	    temp[10] = state.corrRelVort.at(i); temp[11] = state.omeg.at(i); 
	    temp[12] = abs(state.omeg.at(i) - 1);
            temp2[0] = angleQ.at(i); temp2[1] = eigenR.at(i);
	    temp2[2] = state.corrVort_error.at(i);
	    temp2[3] = state.corrRelVort.at(i);

	    if( temp2[2] > pow(10, -3.0)){
	      corr.push_back(temp2);
	    }

	    particles.push_back(temp);
 
	    vecFerror.push_back(errorF);

	    PR_TIMERS("main::generateVec_Pwrite");
           }

	 }

         {
         PWrite<DIM+11>(particles,varnames,testname,k);
         PWrite<DIM+2>(corr, varnames2, filen, k);
         PR_TIMERS("main::pwrite");

         }  

	 {
          errorU.push_back(state.sumU);
          errorVort.push_back(state.sumVort);
          errorVg.push_back(state.sumUg );
          errorU1.push_back(state.sumU1);
          errorVort1.push_back(state.sumVort1);
          errorVg1.push_back(state.sumUg1);
          errorVortm.push_back(state.sumVortm);
	  PR_TIMERS("main::add_error_sum");
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
          rk4.advance(time,dt,state);
	  PR_TIMERS("main::RK4");
	  }
	  {
	  f.update_X(state);
	  PR_TIMERS("main::update_f");
	  }
	  {
	   rk4_f.advance(time, dt, f);
	   PR_TIMERS("main::RK4_f");
	  }

          PR_TIMERS("main::time_loop");
	  

          
    }

  {
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


      for(unsigned int i = 0; i < t_remap.size(); i++){

          f9 << t_remap.at(i) << " " << num_p.at(i)  <<  endl;

      }

   }

   PR_TIMERS("main::write_error_curve");
  }
PR_TIMER_REPORT();

  

  
}
