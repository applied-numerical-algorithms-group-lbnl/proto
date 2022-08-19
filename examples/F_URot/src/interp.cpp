#include "Proto.H"
#include <iostream>
#include <particle.hpp>
#include <vector>
#include "interp.H"
#include <cmath>

using namespace Proto;


void interp::W22( BoxData<double>& A, const vector<particle> p, const double h_to, const double h_from, const int Np){

    Point r;
    int px, py;
    Point i_to;
    Point i_from;
    double x_to[2], x_from[2];
    double W22[2];
    Box Bsten( Point({-2, -2}), Point({2,2}) );
    int n;
    double d[2];
    double strength;
  //  n = L*h_to;
  //  Box B_to( Point::Zeros(), Point( n,n ) );
  //  BoxData<double> A(B_to);
//   cout << h_to << " " << h_from << endl;
    A.setToZero();
    
     for( unsigned int k = 0; k < Np; k++){

      //  cout << p[k].x[0] << " " << p[k].x[1] << endl;

       px = static_cast<int> ( round(p[k].x[0] / h_to) );
       py = static_cast<int> ( round(p[k].x[1] / h_to) );

       i_from = Point({px, py} );

       x_from[0] = p[k].x[0];
       x_from[1] = p[k].x[1];

  
       for( auto i = Bsten.begin(); i != Bsten.end(); i++){


         r = i.operator*();
         i_to = r+i_from;

         x_to[0]=i_to[0]*h_to; x_to[1] = i_to[1]*h_to;


         for(int j =0; j < DIM; j++){

           d[j] = ( abs((x_to[j] - x_from[j]) / h_to ) );

           if((d[j] < 1.0)){

             W22[j] = 1 - pow(d[j], 2.0) - 9.0/2.0*pow( d[j], 3.0) + 15.0/2.0*pow(d[j],4.0) - 3.0*pow( d[j], 5.0);

           } else if( (d[j] >= 1.0) && (d[j] < 2.0) ){

             W22[j] = -4 +18.0*d[j] -29.0*pow(d[j], 2.0) + 43.0/2.0*pow( d[j],3.0) - 15.0/2.0*pow( d[j], 4.0) + pow( d[j] , 5.0);

           } else{

             W22[j] = 0;
           }

         } 

     
          A( i_to ) += p[k].strength *pow(h_from/h_to, 2.0) * ( W22[0]*W22[1]);



      }   

          
   }



}


void interp::W22p(BoxData<double>& A,  vector<double>& p, const vector<particle> X, const double h_to, const double h_from, const int Np){
	

    Point r;
    int px, py;
    Point ig;
    Point i_closest;
    double x_to[2], x_from[2];
    double sum;
    double W22[2];
    double d[2];
    Box Bsten( Point({-2, -2}), Point({2,2}) );
    int n;
    double strength;
    p.clear();    

   for(int k = 0; k < Np; k++){

     x_to[0] = X.at(k).x[0]; x_to[1] = X.at(k).x[1];
   
     px = static_cast<int> ( round(x_to[0] / h_from) );
     py = static_cast<int> ( round(x_to[1] / h_from) );

     i_closest = Point({px, py} );

 
     for( auto i = Bsten.begin(); i != Bsten.end(); i++){
          
           r = i.operator*();
           ig = r+i_closest;
	   
	   strength = *A.data(ig);
	   x_from[0] = ig[0]*h_from; x_from[1] = ig[1]*h_from;
	     
           for(int j =0; j < DIM; j++){

             d[j] = ( abs((x_from[j] - x_to[j]) / h_from ) );

             if((d[j] < 1.0)){

                W22[j] = 1 - pow(d[j], 2.0) - 9.0/2.0*pow( d[j], 3.0) + 15.0/2.0*pow(d[j],4.0) - 3.0*pow( d[j], 5.0);

             } else if( (d[j] >= 1.0) && (d[j] < 2.0) ){

                W22[j] = -4 +18.0*d[j] -29.0*pow(d[j], 2.0) + 43.0/2.0*pow( d[j],3.0) - 15.0/2.0*pow( d[j], 4.0) + pow( d[j] , 5.0);

              } else{

                W22[j] = 0;

              }

            }

                  


            sum += strength * ( W22[0]*W22[1]);


        }


            p.push_back(sum);
	    sum = 0;

      }


   }





