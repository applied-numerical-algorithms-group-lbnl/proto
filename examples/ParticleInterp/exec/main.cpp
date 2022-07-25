#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include <memory>
#include <stdio.h>
#include <fstream>
#include "Proto.H"
#include <iostream>
#include <cstring>
#include <memory>
//#include "LevelMultigrid.H"
//#include "InputParser.H"
#include "particle.hpp"

using namespace std;
using namespace Proto;
//inline 

int main(int argc, char* argv[])
{
    double nx =200;
    double ny = 200;
    int boxSize = 64;
    int L = 1;
    double hg;
    double hp;
    int Np = 0;
    particle m_p;
    vector<particle> p;
    double w;
    Point ip;
    Point ig;
    double xp[2];
    Point r;
    double r2[2];
    int px, py;   
    double tol = pow( 1.0, -6.0);
    double W20_Y, W20_X;
    double sum;
    double xg[2];
   //Define initial grid
    Box Bp(Point::Zeros(), Point(nx-1,ny-1));



   //Define second grid--refinement
    Box Bg(Point::Zeros(),Point(1.2*nx - 1, 1.2*ny -1 ) );
 
    hg = L/(1.2*nx -1);

    BoxData<double> A(Bg);

    BoxData<double> IC(Bp);
    //Define hp
    hp = L/(nx-1);

   IC.setToZero();

    //Initalize particles and generate particle vector
    for( auto iter = Bp.begin(); iter != Bp.end(); iter++)
    {       
            r = iter.operator*();
	    
	    m_p.x[0] = r[0]*hp;
	    m_p.x[1] = r[1]*hp;
	    

	    //Radial profile
            w = 1 - 0.5*(1-erf( (pow(0.3,2) - pow( (m_p.x[0]- 0.5), 2) -  pow((m_p.x[1]-0.5), 2) ) / (0.02) ));
            
//	   Point r2 = r;

	    IC(r)+= w;
	    if ( w > 0 ){

               m_p.strength = w;
	       p.push_back(m_p);
	       Np++;

            }

  }

 
  //Interpolate to refined grid
//  Npp = 0;
   A.setToZero();

   Box Bsten( Point({-1, -1}), Point({1,1}) );

 cout << hg <<" \n";
   for( unsigned int k = 0; k < Np; k++){
     

     px = static_cast<int> ( floor(p[k].x[0] / hg) ); 
     py = static_cast<int> ( floor(p[k].x[1] / hg) );

     ip = Point({px, py} );
  
     //  ip[1] = floor( p[k].x[1] / hg );

     xp[0] = p[k].x[0];
     xp[1] = p[k].x[1];

// cout << xp[0] << " " << xp[1] << "\n";
     //Iterate over spline stencil of contributing points
     for( auto i = Bsten.begin(); i != Bsten.end(); i++){

       r = i.operator*();

       ig = r+ip;

       xg[0]=ig[0]*hg; xg[1] = ig[1]*hg;


    // cout << xg[0] << " " << xg[1] << "\n";
     //  cout << ig << "\n";
     
       W20_X = ( abs((xg[0] - xp[0]) / hg ) );
       W20_Y = ( abs((xg[1] - xp[1]) / hg ) );

     //  cout << W20_X << " " << W20_Y << "\n";
       if((W20_X < 1)){
         W20_X = (1 - abs((xg[0] - xp[0]) / hg ) );

       } 
       else if (W20_X < tol){
         W20_X = 0;
       }
       else{

	W20_X = 0;
       }


       if((W20_Y < 1) ){
         W20_Y = (1 - abs((xg[1] - xp[1]) / hg ) );
       }
       else if (W20_Y < tol) {

	 W20_Y = 0;

       }
       else{
	 W20_Y = 0;
       }


   //   cout << W20_X << " " << W20_Y << "\n";

//       sum +=  p[k].strength * hp/hg * ( W20_X*W20_Y);
       A( ig ) += p[k].strength * hp/hg * ( W20_X*W20_Y);

       

     }

    
   }

double strength;
Point r3;

  IC.setToZero();

 for( auto it = Bg.begin(); it != Bg.end(); it++){
     r3 = it.operator*();
     
   px = static_cast<int> ( floor(r3[0]*hg / hp));
   py = static_cast<int> ( floor(r3[1]*hg / hp) );

   ip = Point({px,py});
     // ip[1] = floor( p[k].x[1] / hg );

     xp[0] = r3[0]*hg;
     xp[1] = r3[1]*hg;
     strength = *A.data(r3);
//    cout << ip << "\n";
    if( strength > 0.000001){

//     cout << "in" << "\n";
	     
	     // cout << xp[0] << " " << xp[1] << "\n";
     //Iterate over spline stencil of contributing points
     for( auto i = Bsten.begin(); i != Bsten.end(); i++){

       r = i.operator*();

       ig = r+ip;

       xg[0]=ig[0]*hp; xg[1] = ig[1]*hp;


    // cout << xg[0] << " " << xg[1] << "\n";
     //  cout << ig << "\n";

       W20_X = ( abs((xg[0] - xp[0]) / hp ) );
       W20_Y = ( abs((xg[1] - xp[1]) / hp ) );

     //  cout << W20_X << " " << W20_Y << "\n";
       if((W20_X < 1) ){
         W20_X = (1 - abs((xg[0] - xp[0]) / hp ) );

       }
       else if (W20_X < tol){
         W20_X = 0;
       }
       else{

        W20_X = 0;
       }

      if((W20_Y < 1) ){
         W20_Y = (1 - abs((xg[1] - xp[1]) / hp ) );
       }
       else if (W20_Y < tol) {

         W20_Y = 0;

       }
       else{
         W20_Y = 0;
       }


//       cout << W20_X << " " << W20_Y << "\n";

//       sum +=  p[k].strength * hp/hg * ( W20_X*W20_Y);
       IC( ig ) += strength * hg/hp * ( W20_X*W20_Y);

     }


     cout << IC(ig) <<"\n";
  }

   }

                       







//    for (int iter = 0; iter < maxIter; iter++)
//    {
//        PR_TIMERS("MG top level");
//        mg.vCycle(phi,rho);



int n = 0;
 #ifdef PR_HDF5
        HDF5Handler h5;
     h5.writePatch(hp,IC, "BacktoIC_Test%i.hdf5", n);
     h5.writePatch(hg,A, "VF_Test%i.hdf5", n);

 #endif
    
//    }

//    PR_TIMER_REPORT();
//#ifdef PR_MPI
//    MPI_Finalize();
//#endif
}

