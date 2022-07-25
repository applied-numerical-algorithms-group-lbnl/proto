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

#define NUMCOMPS DIM+2

using namespace std;
using namespace Proto;


void PrintFile( string fileout, vector<particle> X , double time){

     ofstream f2(fileout);
     string t;
     string p;
     string p2;

     p2 = "# y";
     p = "# Particles";
     t = "# TIME ";
     t.append( to_string( time ) );

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

  int nx = 10;
  int ny = 10;
  Box Bg(Point::Zeros(), Point({nx-1, ny-1}));
  vector<particle> X;
  double dt = 0.1;
  double time = 0;
  double tstop = 1;
  int maxStep;
  Point r;
  double h;
  double L = 1;
  particle p;
  string solnfile;
  string soln_prefix = "SolidRot";
  maxStep = tstop/dt;

  h = L/(nx-1);


  //Initialize Position and Vorticity 
  for( auto it = Bg.begin(); it != Bg.end(); it++){

    r = it.operator*();
    
    p.x[0] = r[0]*h;
    p.x[1] = r[1]*h;
    p.strength = -1;

    X.push_back(p);
    
  }

  State state(X,dt,h);
 
  RK4<State,F,DX> rk4;
  
  

 // Advance RK4 & print to file 
  for (int k = 0;(k < maxStep);k++)
  {


	 
         if (k <10){
             solnfile = soln_prefix +"00";
             }
         else if ( k > 9 && k < 100){
             solnfile = soln_prefix +"0";}
         else{
             solnfile = soln_prefix;
         }
	  solnfile.append( to_string( k ) );
          solnfile.append(".curve");
          PrintFile( solnfile, state.X, time );
          time += dt;
    
          rk4.advance(time,dt,state);
   
//  PR_TIMER_REPORT();
    }
}
