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
  PR_TIME("main");
  int nx = 128;
  int ny = 128;
  Box Bg(Point::Zeros(), Point({nx-1, ny-1}));
  vector<particle> X;
  double dt = 0.0421875;
  double time = 0, timef = 0;
  vector<double> t, angleError, t_remap;
  vector<int> num_p;
  vector<array<double,DIM+13> > particles;
  vector<array<double, DIM+2>> corr;
  deform d;
  double tstop = 20;
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
  double sumVort = 0;
  double sumUg = 0;
  vector<double> u,v, angleQ, eigenR;
  vector<deform> Q, R2, vecFerror;
  deform errorF;
  BoxData<double> vort(Bg);
  BoxData<double> Vortg(Bg);
  Box Bp(Point::Zeros(), Point({4*(nx-1), 4*(nx-1)}));
  BoxData<double, DIM> F_result;
  vector<double> errorU,errorU1,errorVg, errorVg1, errorVort1,errorVort, errorpsi;
  string q,e;
  double Np = 0;
  double wp = 0;

  
  string testname="test";
  string filen = "corr";
  const char *const varnames[] = {"x_1", "x_2", "angleQ", "eigenR", "F11",
                                  "F12", "F21",
                                  "F22", "errorF11",
                                  "errorF12", "errorF21", "errorF22", "vorticity", "absvorterror", "relvorterror"}; 
  const char *const varnames2[] = {"angleQ", "eigenR", "absvorterror", "relvorterror"};

  h = L/(nx-1);
  hp = L/( 4*(nx-1) );
  double h2 = L/nx;
  Ug.setToZero(); Vg.setToZero();
  Vortg.setToZero(); vort.setToZero();


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

   //Initialize Position and Vorticity
  for( auto it = Bp.begin(); it != Bp.end(); it++){

    r = it.operator*();

    p.x[0] = r[0]*hp;
    p.x[1] = r[1]*hp;
 //  cout << p.x[0] << " " << p.x[1] << endl;

    z = ( pow( pow((p.x[0] - c),2.0) + pow((p.x[1] - c),2.0), 0.5 ));


    if( z <= 1){
    
    	    
       p.strength = pow( ( 1 - pow(z, 2.0) ), 7.0 );

       if( p.strength > pow(10, -7.0)){
       Np++;

       wp += p.strength*pow(hp, 2.0);
       X.push_back(p);
       }
    }
  }

num_p.push_back(Np);

#ifdef PR_HDF5
     HDF5Handler h5;
     h5.writePatch(h,Ug,{"U"}, "exact_U%i.hdf5", 0);
     h5.writePatch(h,Vg,{"V"}, "exact_V%i.hdf5", 0);

#endif

  State state(X,dt,hp, h, L, wp);
  deformation f(X); 
  RK4<State,F,DX> rk4;
  RK4<deformation, RHS, DF> rk4_f;
  int step = 0.16875/dt;
  int tpt;

  double sumU1 = 0, sumVort1 = 0, sumUg1 = 0;
  sumU = 0; sumVort = 0; sumUg = 0;
  Vortg.setToZero();

  state.getVelocity(Vortg,errorVg, errorpsi,u,v,Np,0,state);

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

  for(auto it = Bg.begin(); it != Bg.end(); it++){
      r = it.operator*();
      sumVort +=  pow( *Vortg.data(r) - *vort.data(r), 2.0)*pow(h, 2.0);
      sumVort1 +=  abs( *Vortg.data(r) - *vort.data(r))*pow(h, 2.0);



  }

  sumU = pow(sumU, 0.5);
  sumVort = pow(sumVort, 0.5);
  errorU.push_back(sumU); errorVort.push_back(sumVort);
  errorU1.push_back(sumU1); errorVort1.push_back(sumVort1);
  errorVg.push_back(state.sumUg); errorVg1.push_back(state.sumUg1);
 

  t.push_back(0); t_remap.push_back(0.0);
  array<double, DIM+13> temp;
  array<double, DIM+2> temp2;


  for( int k = 0;(k < maxStep);k++)
  {
       //  cout << time << " time " << endl;
         f.QR(Q,R2, angleQ, eigenR);

//	 cout << "got QR" << endl;
         vecFerror.clear(); angleError.clear(); particles.clear(); corr.clear();

       for( int i = 0; i <  state.X.size(); i++){
 
//      F11 = pow(f.F.at(i).F[0][0], 2.0) + pow(f.F.at(i).F[1][0], 2.0);
//      F21 = f.F.at(i).F[0][0]*f.F.at(i).F[0][1] + f.F.at(i).F[1][1]*f.F.at(i).F[1][0];
//      F22 = pow(f.F.at(i).F[1][1], 2.0) + pow(f.F.at(i).F[0][1], 2.0);
            
//   if( i == 0){
//   if( (abs(F11-1) > pow(10, -7.0)) || (abs(F21) > pow(10, -70)) || (abs(F22 - 1) > pow(10, -7.0) ) ){

//     cout << "FTF =/= I" << endl;
//            cout << F11 << " " << F21 << endl;
//            cout << F21 << " " << F22 << endl;	    

	//    }
	    

//   if( abs(f.F.at(i).F[0][0]*f.F.at(i).F[1][1] - f.F.at(i).F[1][0]*f.F.at(i).F[0][1] - 1) > pow(10.0, -7.0)){

//           cout << "det(F) =/= 1 " << abs(F11*F22 - F21*F21-1) << endl;

//   }
//	    }

	    
	    z2 = pow( state.X.at(i).x[0]-1.5, 2.0) + pow(state.X.at(i).x[1]-1.5, 2.0);   

	    x = state.X.at(i).x[0]; y = state.X.at(i).x[1];
            z = ( pow( pow((x - c),2.0) + pow((y - c),2.0), 0.5 ));
    
	    if( z2 > pow( 10.0, -10.0) ){
//          cout << x << " " << y << endl;
	    if( abs(z) > 1){

              dudx = -1*(x-c)*(y-c)/ (8*pow( z2, 2.0)) ; 
	      dudy = -1*(pow(y-c, 2.0) - pow(x-c, 2.0) ) / (16*pow(z2, 2.0) );
	      dvdx = -1*(pow(y-c, 2.0) - pow(x-c, 2.0) ) / (16*pow(z2, 2.0) );
	      dvdy = (x-c)*(y-c)/ (8*pow( z2, 2.0) ); 

	    }else{

              dudx = ((x-c)*(y-c)*(1 - pow(1-z2, 8.0) ))/(8.0*pow(z2, 2.0) ) - ((x-c)*(y-c)*pow( 1 - z2, 7.0))/z2; dudx *= -1.0;
	      dudy = (pow(y-c, 2.0)*(1 - pow( 1-z2, 8.0)))/(8*pow(z2, 2.0)) - (1 - pow( 1-z2, 8.0))/(16*z2) - (pow(y-c, 2.0)*pow(1-z2, 7.0))/z2; dudy *= -1.0;
	      dvdx = (-1*pow(x-c, 2.0)*(1 - pow( 1-z2, 8.0)))/(8*pow(z2, 2.0)) + (1 - pow( 1.0-z2, 8.0))/(16.0*z2) + (pow(x-c, 2.0)*pow(1-z2, 7.0))/z2; dvdx *= -1.0;
	      dvdy = (-1*(x-c)*(y-c)*(1 - pow(1-z2, 8.0) ))/(8.0*pow(z2, 2.0)) + ((x-c)*(y-c)*pow( 1 - z2, 7.0))/z2; dvdy *= -1.0;


	    }

	    if( (abs(dudx) < pow(10.0, -14.0)) && (abs(dvdy) < pow(10, -14.0) )){

	      c_dudy= complex<double> ( dudy, 0.0);
              c_dvdx= complex<double> ( dvdx, 0.0);

	      if( dudy*dvdx > 0 ){
	        lambda1= complex<double>(sqrt( dudy*dvdx), 0.0) ;
	        lambda2= complex<double>(-1.0*sqrt( dudy*dvdx ), 0.0);
	      }else{
                
	        lambda1=complex<double>(0.0, -sqrt(abs(dudy*dvdx))); 
		lambda2= complex<double>(0.0,  sqrt(abs(dudy*dvdx)));


	      }
	     
	      if(dvdx > 0){
	        sqrtv= complex<double> (sqrt(dvdx), 0.0);
	      }else{
	        sqrtv= complex<double> (0, sqrt(abs(dvdx)));	     
	      }

	      if(dudy > 0){
                sqrtu= complex<double> (sqrt(dudy), 0.0);
              }else{
                sqrtu= complex<double> (0, sqrt(abs(dudy)));
              }

	      F11 = (0.5)*( exp(lambda1*(timef)) + exp(lambda2*(timef)) ) ; F12 = (0.5)*sqrtu/sqrtv*(exp(lambda2*(timef)) - exp(lambda1*(timef)));
              F21 = (0.5)*sqrtv/sqrtu*(exp(lambda2*(timef)) - exp(lambda1*(timef))); F22 = (0.5)*( exp(lambda1*(timef)) + exp(lambda2*(timef)) );

	      errorF.F[0][0] = abs( f.F.at(i).F[0][0] - F11.real() );
              errorF.F[0][1] = abs( f.F.at(i).F[0][1] - F12.real() );
              errorF.F[1][0] = abs( f.F.at(i).F[1][0] - F21.real() );
              errorF.F[1][1] = abs( f.F.at(i).F[1][1] - F22.real() );
           // cout << "case 0 " << F11 << " " << F12 << " " << F21 << " " << F22 << endl;
	    } else if( (abs(dudy) < pow(10.0, -14.0)) && (abs(dvdx) < pow(10, -14.0))){

              F11 = exp(-abs(dudx)*timef); F12 = 0;
	      F21 = 0; F22 = exp(abs(dudx)*timef);

	      errorF.F[0][0] = F11.real(); //abs( f.F.at(i).F[0][0] - F11.real() );
              errorF.F[0][1] = F12.real(); //abs( f.F.at(i).F[0][1] - F12.real() );
              errorF.F[1][0] = F21.real(); //abs( f.F.at(i).F[1][0] - F21.real() );
              errorF.F[1][1] = F22.real(); //abs( f.F.at(i).F[1][1] - F22.real() );
//	      cout << "case 1 " << F11 << " " << F12 << " " << F21 << " " << F22 << endl;

	    } else if( (abs(dudy) < pow(10.0, -14.0) ) ){

              F11 = exp(abs(dudx)*timef); F12 = 0;
              F21 = dvdx/(2*dudx)*(exp(abs(dudx)*timef) - exp(-abs(dudx)*timef)); F22 = exp(-abs(dudx)*timef);

              errorF.F[0][0] = F11.real(); //abs( f.F.at(i).F[0][0] - F11.real() );
              errorF.F[0][1] = F12.real(); //abs( f.F.at(i).F[0][1] - F12.real() );
              errorF.F[1][0] = F21.real(); //abs( f.F.at(i).F[1][0] - F21.real() );
              errorF.F[1][1] = F22.real(); //abs( f.F.at(i).F[1][1] - F22.real() );

	      cout << "case 2 " << F11 << " " << F12 << " " << F21 << " " << F22 << endl;


	    }else if( (abs(dvdx) < pow(10.0, -14.0))){


	      F11 = exp(abs(dudx)*timef); F12 =dudy/(2*dudx)*(exp(abs(dudx)*timef) - exp(-abs(dudx)*timef)) ;
              F21 = 0; F22 = exp(-abs(dudx)*timef);

              errorF.F[0][0] = F11.real(); //abs( f.F.at(i).F[0][0] - F11.real() );
              errorF.F[0][1] = F12.real(); //abs( f.F.at(i).F[0][1] - F12.real() );
              errorF.F[1][0] = F21.real(); //abs( f.F.at(i).F[1][0] - F21.real() );
              errorF.F[1][1] = F22.real(); //abs( f.F.at(i).F[1][1] - F22.real() );

	      cout << "case 3 " << F11 << " " << F12 << " " << F21 << " " << F22 << endl;

	    }else{

             c_dudy = complex<double> ( dudy, 0.0);
             c_dvdx = complex<double> ( dvdx, 0.0);
	     c_dudx = complex<double> ( dudx, 0.0);
             c_dvdy = complex<double> ( dvdy, 0.0);


	     if( (pow(dudx, 2.0)+dudy*dvdx) > 0){

               lambda1 = complex<double> (-1.0*( sqrt( pow(dudx, 2.0)+dudy*dvdx)), 0.0 );
	       lambda2 = complex<double> ( sqrt( pow(dudx, 2.0) + dudy*dvdx ),0.0 );
	     }else{

	       lambda1= complex<double> (0.0, -1.0*sqrt(abs(pow(dudx, 2.0)+dudy*dvdx)));
	       lambda2= complex<double> (0.0, sqrt(abs(pow(dudx, 2.0)+dudy*dvdx)));
	     }

	     if( abs(-1.0*dudx*dvdy+dudy*dvdx) < pow(10, -6.0))
		     cout << "small lambda"  << lambda1 << " " << lambda2 << endl;

	       alpha = (c_dudx + lambda1)/dvdx;
	       beta  = (c_dudx + lambda2)/dvdx;
             
	     complex<double> exp2, exp1, ab, a_b;
             double alpha_r, beta_r,F11_r, F12_r, F21_r, F22_r;

	     alpha_r = (dudx+ lambda1.real())/dvdx;
	     beta_r = (dudx + lambda2.real())/dvdx;

	     ab = alpha*beta; a_b = 1.0/(alpha-beta);
	     exp2 = exp(lambda2*timef); exp1 = exp(lambda1*timef); 
	     F11 = a_b*( alpha*exp1 - beta*exp2 );
	     F12 = ab*a_b*( exp2 - exp1 );
	     F21 = a_b*(exp1 - exp2 );
	     F22 = a_b*( alpha*exp2 -  beta*exp1 );
	     
	     errorF.F[0][0] = abs( f.F.at(i).F[0][0] - F11 );
             errorF.F[0][1] = abs( f.F.at(i).F[0][1] - F12 );
             errorF.F[1][0] = abs( f.F.at(i).F[1][0] - F21 );
             errorF.F[1][1] = abs( f.F.at(i).F[1][1] - F22 );

	    } 
	    temp[0] = state.X.at(i).x[0]; temp[1] = state.X.at(i).x[1];
	    temp[2] = angleQ.at(i); temp[3] =eigenR.at(i);
	    temp[4] = f.F.at(i).F[0][0]; temp[5] = f.F.at(i).F[0][1];
            temp[6] = f.F.at(i).F[1][0]; temp[7] = f.F.at(i).F[1][1];
	    temp[8] = errorF.F[0][0]; temp[9] = errorF.F[0][1];
	    temp[10] = errorF.F[1][0]; temp[11] = errorF.F[1][1];
	    temp[12] = state.X.at(i).strength; temp[13] = state.corrVort_error.at(i);
	    temp[14] = state.corrRelVort.at(i);
            temp2[0] = angleQ.at(i); temp2[1] = eigenR.at(i);
	    temp2[2] = state.corrVort_error.at(i);
	    temp2[3] = state.corrRelVort.at(i);

	    if( temp2[2] > pow(10, -3.0)){
	      corr.push_back(temp2);
	    }

	    particles.push_back(temp);
 
	    vecFerror.push_back(errorF);
	 }

	 }

         PWrite<DIM+13>(particles,varnames,testname,k);
         PWrite<DIM+2>(corr, varnames2, filen, k);


  
          errorU.push_back(state.sumU);
          errorVort.push_back(state.sumVort);
          errorVg.push_back(state.sumUg );
          errorU1.push_back(state.sumU1);
          errorVort1.push_back(state.sumVort1);
          errorVg1.push_back(state.sumUg1);
 

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


          rk4.advance(time,dt,state);
	  f.update_X(state);
	//  cout << "F update" << endl;
	  rk4_f.advance(time, dt, f);
	//  cout << "F advance" << endl;

        //f( (k+1)%5 == 0){
        //  state.remap();
    	//    f.remap( state);
    	//    timef = 0;
        // }
	  

          
    }

  string errorfile = "errorVelocity.curve";
  string errorfile2 = "errorVorticity.curve";
  string errorfile3 = "errorVg.curve";
  ofstream f3(errorfile);
  ofstream f4(errorfile2);
  ofstream f5(errorfile3);
  string errorfile4 = "errorVelocity1.curve";
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


  e = "# Velocity1";


  if( f6.is_open() ) {

      //format in scientific notation
       f6.setf(std::ios::scientific, std::ios::floatfield);
       f6.precision(4);


      f6 << q << endl;


      f6 << e << endl;


      for(unsigned int i = 0; i < t.size(); i++){

          f6 << t.at(i) << " " << errorU1.at(i)  <<  endl;

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

PR_TIMER_REPORT();

  

  
}
