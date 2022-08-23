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


#define NUMCOMPS DIM+2

using namespace std;
using namespace Proto;


void PrintFile( string fileout,string fileout2, vector<deform> f, vector<particle> X, vector<deform> errorf, vector<double> angleQ, vector<double> angleError, vector<double> eigenR, double time){
     
 //    ofstream f1(fileout2);
     ofstream f2(fileout);
     string t;
     string p;
     string p2;
     string a,d, q, r;

     d = "#F"; a = "#angle error";
     p2 = "#error in F ";
     r = "#max eigenvalue of R";
     q = "#angle of Q";
     p = "# Particles";
     t = "# TIME ";
     t.append( to_string( time ) );

//   if(f1.is_open() ){

 //       f1.setf(std::ios::scientific, std::ios::floatfield);
 //       f1.precision(7);

// f1 << t << endl;
// f1 << d << endl;

 //        for(unsigned int i = 0; i < f.size(); i++){
//
 //          f1 << X.at(i).x[0] << " " << X.at(i).x[1] << endl;
 //          f1 << f.at(i).F[0][0] << " " << f.at(i).F[0][1] << endl;
//    f1 << f.at(i).F[1][0] << " " << f.at(i).F[1][1] << endl;
//    f1 << endl;

//         }
//  f1 << endl;
//  f1 << p2 << endl;
           
//  for(unsigned int i = 0; i < f.size(); i++){
//    f1 << X.at(i).x[0] << " " << X.at(i).x[1] << endl;
//           f1 << errorf.at(i).F[0][0] << " " << errorf.at(i).F[0][1] << endl;
 //          f1 << errorf.at(i).F[1][0] << " " << errorf.at(i).F[1][1] << endl;
//           f1 << endl;

//   }

//  f1 << endl;
//  f1 << q << endl;
 //        for(unsigned int i = 0; i < f.size(); i++){
 //          f1 << X.at(i).x[0] << " " << X.at(i).x[1] << endl;
 //          f1 <<  angleQ[i] << endl;
//    f1 << endl;

//          }

//         f1 << endl;
//         f1 << r << endl;
//         for(unsigned int i = 0; i < f.size(); i++){
//           f1 << X.at(i).x[0] << " " << X.at(i).x[1] << endl;
//           f1 <<  eigenR[i] << endl;
//           f1 << endl;

//          }

//  f1 << endl;
//         f1 << a << endl;
//         for(unsigned int i = 0; i < f.size(); i++){
//           f1 << X.at(i).x[0] << " " << X.at(i).x[1] << endl;
//           f1 <<  angleError.at(i) << endl;
//           f1 << endl;

//          }





//   }


     if(f2.is_open()){
       f2.setf(std::ios::scientific, std::ios::floatfield);
           f2.precision(4);
   
          //write solution to file
     
           f2 << t << endl;
            
                
           f2 << p << endl;
               

           for(unsigned int i = 0; i < X.size(); i++){

		   f2 << X.at(i).x[0] << " " << X.at(i).x[1] <<  endl;

            }


     }

}



/***/
int main(int argc, char* argv[])
{

  int nx = 32;
  int ny = 32;
  Box Bg(Point::Zeros(), Point({nx-1, ny-1}));
  vector<particle> X;
  double dt = 0.16875;
  double time = 0;
  vector<double> t, angleError;
  vector<array<double,DIM+11> > particles;
  deform d;
  double tstop = 20;
  double maxStep;
  Point r;
  double h;
  double R = 1;
  double hp;
  double F11, F21, F22;
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
  double z, z2;
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
  vector<double> errorU,errorVg, errorVort, errorpsi;
  string q,e;
  double Np = 0;
  double wp = 0;

  string testname="test";
  const char *const varnames[] = {"x_1", "x_2", "angleQ", "eigenR", "F11",
                                  "F12", "F21",
                                  "F22", "errorF11",
                                  "errorF12", "errorF21", "errorF22", "vorticity"};

  h = L/(nx-1);
  hp = L/( 4*(nx-1) );
  double h2 = L/nx;
  Ug.setToZero(); Vg.setToZero();
  Vortg.setToZero(); vort.setToZero();


 for( auto it = Bg.begin(); it != Bg.end(); it++){

    r = it.operator*();

    p.x[0] = r[0]*h;
    p.x[1] = r[1]*h;
 //  cout << p.x[0] << " " << p.x[1] << endl;

    z = ( pow( pow((p.x[0] - c),2.0) + pow((p.x[1] - c),2.0), 0.5 ));

    if( z <= R){
        vort(r) += pow( ( 1 - pow(z, 2.0) ), 7.0 );

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
       Np++;

       wp += p.strength*pow(hp, 2.0);
       X.push_back(p);
    }
  }



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


  sumU = 0; sumVort = 0; sumUg = 0;
   Vortg.setToZero();

   state.getVelocity(Vortg,errorVg, errorpsi,u,v,Np,0,state);

 for(int i = 0; i < Np; i++){

       sumU +=  (abs(u[i] - 0.5*(state.X.at(i).x[1] - 1.5)) +abs(v[i] + 0.5*(state.X.at(i).x[0] - 1.5)))*pow(hp,2.0);
    }

    for(auto it = Bg.begin(); it != Bg.end(); it++){
      r = it.operator*();

      sumVort +=  abs( *Vortg.data(r) - *vort.data(r))*pow(h, 2.0);
     }

     errorU.push_back(sumU); errorVort.push_back(sumVort);


     t.push_back(0);
    array<double, DIM+11> temp;



  for( int k = 0;(k < maxStep);k++)
  {

         f.QR(Q,R2, angleQ, eigenR);
         vecFerror.clear(); angleError.clear(); particles.clear();
       //  cout << time << endl;
       for( int i = 0; i <  state.X.size(); i++){
 
   //       F11 = pow(f.F.at(i).F[0][0], 2.0) + pow(f.F.at(i).F[1][0], 2.0);
//   F21 = f.F.at(i).F[0][0]*f.F.at(i).F[0][1] + f.F.at(i).F[1][1]*f.F.at(i).F[1][0];
  //        F22 = pow(f.F.at(i).F[1][1], 2.0) + pow(f.F.at(i).F[0][1], 2.0);
            
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

	    errorF.F[0][0] = abs( f.F.at(i).F[0][0] - cos(time/2) );
	    errorF.F[0][1] = abs( f.F.at(i).F[0][1] - sin(time/2) );
	    errorF.F[1][0] = abs( f.F.at(i).F[1][0] + sin(time/2) );
	    errorF.F[1][1] = abs( f.F.at(i).F[1][1] - cos(time/2) );
                 
	    angleError.push_back( abs(time/2-angleQ.at(i)) );
	 
	    temp[0] = state.X.at(i).x[0]; temp[1] = state.X.at(i).x[1];
	    temp[2] = angleQ.at(i); temp[3] =eigenR.at(i);
	    temp[4] = f.F.at(i).F[0][0]; temp[5] = f.F.at(i).F[0][1];
            temp[6] = f.F.at(i).F[1][0]; temp[7] = f.F.at(i).F[1][1];
	    temp[8] = errorF.F[0][0]; temp[9] = errorF.F[0][1];
	    temp[10] = errorF.F[1][0]; temp[11] = errorF.F[1][1];
	  //  cout << "new temp" << endl;
	    temp[12] = state.X.at(i).strength;

	    particles.push_back(temp);
 
	 //   cout << "added new temp" << endl;
	    vecFerror.push_back(errorF);
	 }

         PWrite<DIM+11>(particles,varnames,testname,k);

//	 cout << " wrote new temp" << endl;

//	 cout << "pwrite" << endl;
 //        if (k <10){
 //              solnfile = soln_prefix +"000";
 //              solnfilef = solnf_prefix +"000";
 //              }
 //          else if ( k > 9 && k < 100){
 //              solnfile = soln_prefix +"00";
 //              solnfilef = solnf_prefix+"00";
 //          }
 //          else if ( k > 99 && k < 1000){
 //             solnfile = soln_prefix +"0";
 //            solnfilef = solnf_prefix+"0";
 //          }
 //          else{
 //              solnfile = soln_prefix;
 //              solnfilef = solnf_prefix;
 //          }
 //           solnfile.append( to_string( k ) );
 //           solnfile.append(".curve");
 //           solnfilef.append( to_string( k) );
 //           solnfilef.append(".txt");
 //           PrintFile( solnfile, solnfilef, f.F, state.X,vecFerror, angleQ, angleError, eigenR, time );
  
          errorU.push_back(state.sumU);
          errorVort.push_back(state.sumVort);
          errorVg.push_back(state.sumUg );
        //if(k%step == 0){
           
            t.push_back(time);

        //}

          time += dt;
	  
//        if( (k+1)%3 == 0){
 
//          state.remap();
            
// }

          rk4.advance(time,dt,state);
	  f.update_X(state);
	  rk4_f.advance(time, dt, f);

         if( (k+1)%5 == 0){

            state.remap();
            f.remap( state);
          }
	  

          
    }


  string errorfile = "errorVelocity.curve";
  string errorfile2 = "errorVorticity.curve";
  string errorfile3 = "errorVg.curve";
  ofstream f3(errorfile);
  ofstream f4(errorfile2);
  ofstream f5(errorfile3);
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


  e = "# Vg";
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
  

  
}
