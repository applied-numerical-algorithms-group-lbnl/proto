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
#include "Proto_WriteBoxData.H"
#include "Proto_Timer.H"
#include "deformation.H"
#include "deform.hpp"

#define NUMCOMPS DIM+2

using namespace std;
using namespace Proto;


void PrintFile( string fileout, string fileout2, vector<deform> f, vector<particle> X , double time){

     ofstream f2(fileout);
     ofstream f1(fileout2);
     string t;
     string p;
     string p2;

     p2 = "# y";
     p = "# Particles";
     t = "# TIME ";
     t.append( to_string( time ) );

     if(f1.is_open() ){

          f1.setf(std::ios::scientific, std::ios::floatfield);
          f2.precision(7);

	   for(unsigned int i = 0; i < f.size(); i++){

             f1 << f.F[0][0] << " " << f.F[0][1] << " " << f.F[1][0] << " " << f.F[1][1] << endl;

           }

     }

     if( f2.is_open() ) {
     
         //format in scientific notation
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
  int ny =32;
  Box Bg(Point::Zeros(), Point({nx-1, ny-1}));
  vector<particle> X;
  double dt = 0.9018228;
  double time = 0;
  vector<double> t;
  double tstop = 20;
  double R = 1;
  int maxStep;
  Point r;
  double h;
  double hp;
  double L = 3;
  BoxData<double> Ug(Bg);
  BoxData<double> Vg(Bg);
  BoxData<double> Ui(Bg);
  BoxData<double> Vi(Bg);
  particle p;
  string solnfile, solnfilef;
  string soln_prefix = "Vort2D";
  string solnf_prefix = "F2D";
  maxStep = tstop/dt;
  double z,z2;
  double sumU = 0;
  double sumVort = 0;
  double sumUg = 0;
  vector<double> u,v;
  BoxData<double> vort(Bg);
  BoxData<double> Vortg(Bg);
  Box Bp(Point::Zeros(), Point({4*(nx-1), 4*(nx-1)}));
  BoxData<double, DIM> F_result;
  vector<double> errorU,errorVg, errorVort, errorpsi;
  string q,e;
  double Np = 0;
  double wp = 0;
  double c = 1.5;
  h = L/(nx-1);
  hp = L/( 4*(nx-1) );
  double h2 = L/nx;
  vector<deform> Q, R2;
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

  


   sumU = 0; sumVort = 0; sumUg = 0;
   Vortg.setToZero();

   state.getVelocity(Vortg,errorVg, errorpsi,u,v,Np,0,state);

#ifdef PR_HDF5
     h5.writePatch(h, vort, "vort%i.hdf5", 0);     
     h5.writePatch(h,Vortg, "vortg%i.hdf5", 0);


#endif


   for(int i = 0; i < Np; i++){

    z =( pow( pow((state.X.at(i).x[0] - c),2.0) + pow((state.X.at(i).x[1] - c),2.0), 0.5 ));

    if( z <= 1){

       sumU +=  (abs(u[i] +1.0/(16.0*pow(z, 2.0) )*(1 - pow( (1 - pow(z, 2.0) ), 8.0) )*(state.X.at(i).x[1] - c)) +abs(v[i] -1.0/(16.0*pow(z, 2.0) )*(1 - pow( (1 - pow(z, 2.0) ), 8.0) )*(state.X.at(i).x[0] - c)))*pow(hp, 2.0);
    }else{
       
       z2 = pow((p.x[0] - c),2.0) + pow((p.x[1] - c),2.0);
        
       sumU += ( abs( u[i] + 1.0/(16.0*z2)*(state.X.at(i).x[1] - c)) + abs( v[i] - 1.0/(16.0*z2) *(p.x[0] - c) ) )*pow(hp, 2.0);

    }
    }

   for(auto it = Bg.begin(); it != Bg.end(); it++){
      r = it.operator*();

      sumVort +=  pow( *Vortg.data(r) - *vort.data(r), 2.0)*pow(h, 2.0);
   
    }

   sumVort = pow(sumVort, 0.5);
     errorU.push_back(sumU); errorVort.push_back(sumVort);
     t.push_back(0);
     


  for (int k = 0;(k < maxStep);k++)
  {
   
          if (k <10){
               solnfile = soln_prefix +"000";
	       solnfilef = solnf_prefix +"000";
               }
           else if ( k > 9 && k < 100){
               solnfile = soln_prefix +"00";
	       solnfilef = solnf_prefix+"00";
	   }
           else if ( k > 99 && k < 1000){
              solnfile = soln_prefix +"0";
	      solnfilef = solnf_prefix+"0";
           }
           else{
               solnfile = soln_prefix;
	       solnfilef = solnf_prefix;
           }
            solnfile.append( to_string( k ) );
            solnfile.append(".curve");
	    solnfilef.append( to_string( k) );
	    solnfilef.append(".txt")
            PrintFile( solnfile, solnfilef, f.F, state.X, time );
   
	  time += dt;
	  
          f.update_X(state);
	
          rk4.advance(time,dt,state);

          rk4_f.advance(time, dt, f); 	 
          f.QR(Q, R2);	 
	
	  if( (k != 0)){
	    errorU.push_back(state.sumU); 
	    errorVort.push_back(state.sumVort);
	    errorVg.push_back(state.sumUg );
            t.push_back(time);
	  }

//	  if( k%step == 0){

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
