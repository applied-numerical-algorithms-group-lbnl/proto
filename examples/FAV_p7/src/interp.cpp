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
    int n = 4.0/h_to;
    Box Bg( Point::Zeros(), Point({n,n}));
    double d[2];
    double strength;
    double c = 0.75;
    double w_to = 0;
    double w_from = 0;
    A.setToZero();
    
     for( unsigned int k = 0; k < Np; k++){

       px = static_cast<int> ( round((p[k].x[0]+c) / h_to) );
       py = static_cast<int> ( round((p[k].x[1]+c) / h_to) );

       i_from = Point({px, py} );

       x_from[0] = p[k].x[0];
       x_from[1] = p[k].x[1];

  
       for( auto i = Bsten.begin(); i != Bsten.end(); i++){


         r = i.operator*();
         i_to = r+i_from;
         x_to[0]=i_to[0]*h_to-c; x_to[1] = i_to[1]*h_to-c;


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

//   for(int i = 0; i < Np; i++){

//     w_from += p[i].strength*pow(h_from, 2.0);

//   }

//   for(auto it = Bg.begin(); it != Bg.end(); it++){

//    r = *it;

//    w_to += *A.data(r)*pow(h_to, 2.0);

//   }

 //    if( abs(w_from-w_to) > pow(10, -16.0) ){

 //      cout << "Issue with conservation in interp from particles to grid " << abs(w_from-w_to) << endl;


 //    }



}


void interp::W44( BoxData<double>& A, const vector<particle> p, const double h_to, const double h_from, const int Np){

    Point r;
    int px, py;
    Point i_to;
    Point i_from;
    double x_to[2], x_from[2];
    double W44[2];
    Box Bsten( Point({-3, -3}), Point({3,3}) );
    double d[2];
    double c = 0.75;
    double strength;
    double w_to = 0;
    double w_from = 0;
    A.setToZero();
    int n = 10;
    double a[n], b[n], g[n];
    double ratio = pow(h_from/h_to, 2.0);

    //Coefficients for interpolation
    a[0] = 1.0; a[1] = 0; a[2] = -5.0/4.0;
    a[3] = 0.0; a[4] = 1.0/4.0; a[5] = -100.0/3.0;
    a[6] = 455.0/4.0; a[7] = -295.0/2.0;
    a[8] = 345.0/4.0; a[9]= -115.0/6.0;

    b[0] = -199.0; b[1] = 5485.0/4.0; b[2] = -32975.0/8.0;
    b[3] = 28425.0/4.0; b[4] = -61953.0/8.0;
    b[5] = 33175/6.0; b[6] = -20685.0/8.0;
    b[7] = 3055.0/4.0; b[8] = -1035.0/8.0;
    b[9] = 115.0/12.0;

    g[0] = 5913.0; g[1] = -89235.0/4.0;
    g[2] = 297585.0/8.0; g[3] = -143895.0/4.0;
    g[4] = 177871.0/8.0; g[5] = -54641.0/6.0;
    g[6] = 19775.0/8.0; g[7] = -1715.0/4.0;
    g[8]  = 345.0/8.0; g[9] = -23.0/12.0;


     for( unsigned int k = 0; k < Np; k++){

       //Find closest point
       px = static_cast<int> ( round((p[k].x[0]+c) / h_to) );
       py = static_cast<int> ( round((p[k].x[1]+c) / h_to) );

       i_from = Point({px, py} );

       x_from[0] = p[k].x[0];
       x_from[1] = p[k].x[1];


       for( auto i = Bsten.begin(); i != Bsten.end(); i++){

       
       	 {
	 //Calculate position of stencil point	 
         r = i.operator*();
         i_to = r+i_from;

         x_to[0]=i_to[0]*h_to-c; x_to[1] = i_to[1]*h_to-c;
	 }

         for(int j =0; j < DIM; j++){

	   //difference between particle location and stencil location	 
           d[j] = ( abs((x_to[j] - x_from[j]) / h_to ) );

	   W44[j] = 0.0;

	   //Generate W44 based on distance d
           if((d[j] < 1.0)){

	     for(int p = (n-1); p > -1; p--){
                 W44[j] = W44[j]*d[j] + a[p];
	     }

           } else if( (d[j] >= 1.0) && (d[j] < 2.0) ){


            for(int p = (n-1); p > -1; p--){
        	  W44[j] = W44[j]*d[j] + b[p];
              }
	   } else if( (d[j] >= 2.0) && (d[j] < 3.0) ){


	       for(int p = (n-1); p > -1; p--){
                  W44[j] = W44[j]*d[j] + g[p];
               }

           } else{

             W44[j] = 0.0;
           }

         }

          {
          //Contribution to grid location 
          A( i_to ) += p[k].strength *ratio * ( W44[0]*W44[1]);
          }


      }


   }

//   for(int i = 0; i < Np; i++){

//     w_from += p[i].strength*pow(h_from, 2.0);

//   }

//   for(auto it = Bg.begin(); it != Bg.end(); it++){

//    r = *it;

//    w_to += *A.data(r)*pow(h_to, 2.0);

//   }

 //    if( abs(w_from-w_to) > pow(10, -16.0) ){

 //      cout << "Issue with conservation in interp from particles to grid " << abs(w_from-w_to) << endl;


 //    }



}



void interp::W22p(BoxData<double>& A,  vector<double>& p, const vector<particle> X, const double h_to, const double h_from, const int Np){
	

    Point r;
    int px, py;
    Point ig;
    Point i_closest;
    double x_to[2], x_from[2];
    double sum, c = 0.75;
    double W22[2];
    double d[2];
    Box Bsten( Point({-2, -2}), Point({2,2}) );
    int n = 4.0/h_from;
    double strength;
    double w_to = 0;
    double w_from = 0;
    Box Bg( Point({2,2}), Point({n-2,n-2}));
    p.clear();    

   for(int k = 0; k < Np; k++){

     x_to[0] = X.at(k).x[0]; x_to[1] = X.at(k).x[1];
   
     px = static_cast<int> ( round((x_to[0]+c) / h_from) );
     py = static_cast<int> ( round((x_to[1]+c) / h_from) );

     i_closest = Point({px, py} );
 
     for( auto i = Bsten.begin(); i != Bsten.end(); i++){
          
           r = i.operator*();
           ig = r+i_closest;
	   
	   strength = *A.data(ig);
	   x_from[0] = ig[0]*h_from-c; x_from[1] = ig[1]*h_from-c;
	     
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


//   for(int i = 0; i < Np; i++){

//     w_to += p[i]*pow(h_to, 2.0);

//   }

//   for(auto it = Bg.begin(); it != Bg.end(); it++){

//    r = *it;

//    w_from += *A.data(r)*pow(h_from, 2.0);

//   }

//     if( abs(w_from-w_to) > pow(10, -16.0) ){

//       cout << "Issue with conservation in interp from grid to particles " << abs(w_from-w_to) << endl;


//     }



}


void interp::W44p(BoxData<double>& A,  vector<double>& p, const vector<particle> X, const double h_to, const double h_from, const int Np){


    Point r;
    int px, py;
    Point ig;
    Point i_closest;
    double x_to[2], x_from[2];
    double sum;
    double W44[2];
    double d[2];
    Box Bsten( Point({-3, -3}), Point({3,3}) );
    double strength;
    double w_to = 0;
    double c = 0.75;
    double w_from = 0;
    p.clear();
    int n = 10;
    double a[n], b[n], g[n];

    a[0] = 1.0; a[1] = 0; a[2] = -5.0/4.0; 
    a[3] = 0.0; a[4] = 1.0/4.0; a[5] = -100.0/3.0;
    a[6] = 455.0/4.0; a[7] = -295.0/2.0;
    a[8] = 345.0/4.0; a[9]= -115.0/6.0;

    b[0] = -199.0; b[1] = 5485.0/4.0; b[2] = -32975.0/8.0;
    b[3] = 28425.0/4.0; b[4] = -61953.0/8.0;
    b[5] = 33175/6.0; b[6] = -20685.0/8.0;
    b[7] = 3055.0/4.0; b[8] = -1035.0/8.0;
    b[9] = 115.0/12.0;

    g[0] = 5913.0; g[1] = -89235.0/4.0;
    g[2] = 297585.0/8.0; g[3] = -143895.0/4.0;
    g[4] = 177871.0/8.0; g[5] = -54641.0/6.0;
    g[6] = 19775.0/8.0; g[7] = -1715.0/4.0;
    g[8]  = 345.0/8.0; g[9] = -23.0/12.0;

   for(int k = 0; k < Np; k++){

     //Calculate grid index of each particle
     x_to[0] = X.at(k).x[0]; x_to[1] = X.at(k).x[1];

     px = static_cast<int> ( round((x_to[0]+c) / h_from) );
     py = static_cast<int> ( round((x_to[1]+c) / h_from) );

     i_closest = Point({px, py} );


     for( auto i = Bsten.begin(); i != Bsten.end(); i++){

	   //Obtain point describing the stencil location relative to particle location
           r = i.operator*();
           ig = r+i_closest;

           strength = *A.data(ig);
           x_from[0] = ig[0]*h_from-c; x_from[1] = ig[1]*h_from-c;

	   //Determine W44 expression using Horner's Rule
           for(int j =0; j < DIM; j++){
    
    	      d[j] = ( abs((x_from[j] - x_to[j]) / h_from ) );

             W44[j] = 0;
             if((d[j] < 1.0)){

	       for(int p = (n-1); p > -1; p--){
	         W44[j] = W44[j]*d[j] + a[p];
               }
               	       
             } else if( (d[j] >= 1.0) && (d[j] < 2.0) ){


               for(int p = (n-1); p > -1; p--){
                 W44[j] = W44[j]*d[j] + b[p];
               }

             } else if( (d[j] >= 2.0) && (d[j] < 3.0) ){


              for(int p = (n-1); p > -1; p--){
                 W44[j] = W44[j]*d[j] + g[p];
              }


             } else{

               W44[j] = 0;
             }

            }


	    //Interpolation contribution
            sum += strength * ( W44[0]*W44[1]);


        }


            p.push_back(sum);
            sum = 0;

      }


//   for(int i = 0; i < Np; i++){

//     w_to += p[i]*pow(h_to, 2.0);

//   }

//   for(auto it = Bg.begin(); it != Bg.end(); it++){

//    r = *it;

//    w_from += *A.data(r)*pow(h_from, 2.0);

//   }

//     if( abs(w_from-w_to) > pow(10, -16.0) ){

//       cout << "Issue with conservation in interp from grid to particles " << abs(w_from-w_to) << endl;


//     }



}



