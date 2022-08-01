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
    double nx = 64;
    double ny = 64;
    int boxSize = 64;
    int L = 1;
    double hg;
    double hp;
    int Np = 0;
    particle m_p;
    vector<particle> p;
    double w;
    double d[2];
    Point ip;
    Point ig;
    double xp[2];
    Point r;
    double r2[2];
    int px, py;   
    double tol = pow( 1.0, -6.0);
    double W22[2];
    double sum;
    double xg[2];
   //Define initial grid
    Box Bp(Point::Zeros(), Point(nx-1,ny-1));



   //Define second grid--refinement
    Box Bg(Point::Zeros(),Point(0.5*(nx-1),0.5*(ny -1)) );
 
    hg = L/(0.5*(nx-1) );

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
            w = 1 - 0.5*(1-erf( (pow(0.25,2.0)- pow( (m_p.x[0]- 0.5), 2.0) -  pow((m_p.x[1]-0.5), 2.0) ) / (0.015) ));
            
//	   Point r2 = r;

	    IC(r)+= w;
	    if ( w > 0 ){

               m_p.strength = w;
	       p.push_back(m_p);
	       Np++;

            }

  }

 
  //Interpolate to refined grid

   A.setToZero();

   Box Bsten( Point({-2, -2}), Point({2,2}) );

   for( unsigned int k = 0; k < Np; k++){
     
 //    cout << p[k].x[0] / hg << " " << p[k].x[1] / hg << endl;
     px = static_cast<int> ( round(p[k].x[0] / hg) ); 
     py = static_cast<int> ( round(p[k].x[1] / hg) );

    
     ip = Point({px, py} );
 //   cout << ip << endl;
     xp[0] = p[k].x[0];
     xp[1] = p[k].x[1];

     //Iterate over spline stencil of contributing points
     for( auto i = Bsten.begin(); i != Bsten.end(); i++){

       r = i.operator*();

       ig = r+ip;

       xg[0]=ig[0]*hg; xg[1] = ig[1]*hg;


        for(int j =0; j < DIM; j++){
     
           d[j] = ( abs((xg[j] - xp[j]) / hg ) );
   
           if((d[j] < 1.0)){

             W22[j] = 1 - pow(d[j], 2.0) - 9.0/2.0*pow( d[j], 3.0) + 15.0/2.0*pow(d[j],4.0) - 3.0*pow( d[j], 5.0);
 
           } else if( (d[j] >= 1.0) && (d[j] < 2.0) ){

             W22[j] = -4 +18.0*d[j] -29.0*pow(d[j], 2.0) + 43.0/2.0*pow( d[j],3.0) - 15.0/2.0*pow( d[j], 4.0) + pow( d[j] , 5.0);
       
       	   } else{

             W22[j] = 0;
	 
	   }

       }  

     if((W22[0] == 1) && (W22[1] == 1) && (p[k].strength > 0.5) ){
      cout << ig << endl;
      cout << d[0] << " " << d[1] << endl;
      cout << W22[0] << " " << W22[1] << endl;
      cout <<  p[k].strength*( W22[0]*W22[1]) << endl;
     }
        A( ig ) += p[k].strength *pow(hp/hg, 1.0) * ( W22[0]*W22[1]);
	

    }

      
 }

cout << hp/hg << endl;

double strength;
Point r3;

IC.setToZero();

 for( int k = 0; k < Np; k++){
   
     
   px = static_cast<int> ( round(p[k].x[0] / hp));
   py = static_cast<int> ( round(p[k].x[1] / hp) );

   ip = Point({px,py});
     
   xp[0] =  p[k].x[0]; xp[1] = p[k].x[1];
  
     //Iterate over spline stencil of contributing points
     for( auto i = Bsten.begin(); i != Bsten.end(); i++){

       r = i.operator*();

       px = static_cast<int> ( round(p[k].x[0] / hg));
       py = static_cast<int> ( round(p[k].x[1] / hg) );


       r3 = Point({px,py});

       ig = r3 + r;

       xg[0]=ig[0]*hg; xg[1] = ig[1]*hg;

       strength = *A.data(ig);

       for(int j =0; j < DIM; j++){

         d[j] = ( abs((xg[j] - xp[j]) / hp ) );

         if((d[j] < 1.0)){

           W22[j] = 1 - pow(d[j], 2.0) - 9.0/2.0*pow( d[j], 3.0) + 15.0/2.0*pow(d[j],4.0) - 3.0*pow( d[j], 5.0);

         } else if( (d[j] >= 1.0) && (d[j] < 2.0) ){

           W22[j] = -4 +18.0*d[j] -29.0*pow(d[j], 2.0) + 43.0/2.0*pow( d[j],3.0) - 15.0/2.0*pow( d[j], 4.0) + pow( d[j] , 5.0);

         } else{

           W22[j] = 0;
         }
 
       }

       IC( ip ) += strength *pow( hg/hp, 1) * ( W22[0]*W22[1]);

      }

 }

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

