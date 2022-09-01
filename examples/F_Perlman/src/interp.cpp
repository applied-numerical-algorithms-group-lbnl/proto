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
    int n = 3/h_to;
    Box Bg( Point::Zeros(), Point({n,n}));
    double d[2];
    double strength;
    double w_to = 0;
    double w_from = 0;
    A.setToZero();
    
     for( unsigned int k = 0; k < Np; k++){


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
    int n = 3/h_to;
    Box Bg( Point::Zeros(), Point({n,n}));
    double d[2];
    double strength;
    double w_to = 0;
    double w_from = 0;
    A.setToZero();

     for( unsigned int k = 0; k < Np; k++){


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
             W44[j] = 1.0 - 5.0/4.0*pow(d[j], 2.0) +1.0/4.0*pow( d[j], 4.0) - 100.0/3.0*pow(d[j], 5.0)+ 455.0/4.0*pow(d[j],6.0) - 295.0/2.0*pow( d[j], 7.0) + 345.0/4.0*pow(d[j], 8.0) - 115.0/6.0*pow(d[j], 9.0);

           } else if( (d[j] >= 1.0) && (d[j] < 2.0) ){

             W44[j] = -199.0 + 5485.0/4.0*d[j] -32975.0/8.0*pow(d[j], 2.0) + 28425.0/4.0*pow( d[j],3.0) - 61953.0/8.0*pow( d[j], 4.0) + 33175.0/6.0*pow( d[j] , 5.0) - 20685.0/8.0*pow(d[j], 6.0) + 3055.0/4.0*pow(d[j], 7.0) - 1035.0/8.0*pow(d[j], 8.0) + 115.0/12.0*pow(d[j], 9.0);

	   } else if( (d[j] >= 2.0) && (d[j] < 3.0) ){

	     W44[j] = 5913.0 - 89235.0/4.0*d[j] + 297585.0/8.0*pow(d[j], 2.0) - 143895.0/4.0*pow( d[j], 3.0) + 177871.0/8.0*pow( d[j], 4.0) - 54641.0/6.0*pow(d[j], 5.0) + 19775.0/8.0*pow(d[j], 6.0) - 1715.0/4.0*pow(d[j], 7.0) + 345.0/8.0*pow(d[j], 8.0) - 23.0/12.0*pow(d[j], 9.0);

           } else{

             W44[j] = 0.0;
           }

         }


          A( i_to ) += p[k].strength *pow(h_from/h_to, 2.0) * ( W44[0]*W44[1]);



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
    double sum;
    double W22[2];
    double d[2];
    Box Bsten( Point({-2, -2}), Point({2,2}) );
    int n = 3/h_from;
    double strength;
    double w_to = 0;
    double w_from = 0;
    Box Bg( Point({2,2}), Point({n-2,n-2}));
    p.clear();    

   for(int k = 0; k < Np; k++){

     x_to[0] = X.at(k).x[0]; x_to[1] = X.at(k).x[1];
   
     px = static_cast<int> ( round(x_to[0] / h_from) );
     py = static_cast<int> ( round(x_to[1] / h_from) );

     i_closest = Point({px, py} );
//cout << k << " " << x_to[0] << " " << x_to[1] << endl;
 
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
    int n = 3/h_from;
    double strength;
    double w_to = 0;
    double w_from = 0;
    Box Bg( Point({2,2}), Point({n-2,n-2}));
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
               W44[j] = 1.0 - 5.0/4.0*pow(d[j], 2.0) +1.0/4.0*pow( d[j], 4.0) - 100.0/3.0*pow(d[j], 5.0)+ 455.0/4.0*pow(d[j],6.0) - 295.0/2.0*pow( d[j], 7.0) + 345.0/4.0*pow(d[j], 8.0) - 115.0/6.0*pow(d[j], 9.0);

             } else if( (d[j] >= 1.0) && (d[j] < 2.0) ){

               W44[j] = -199.0 + 5485.0/4.0*d[j] -32975.0/8.0*pow(d[j], 2.0) + 28425.0/4.0*pow( d[j],3.0) - 61953.0/8.0*pow( d[j], 4.0) + 33175.0/6.0*pow( d[j] , 5.0) - 20685.0/8.0*pow(d[j], 6.0) + 3055.0/4.0*pow(d[j], 7.0) - 1035.0/8.0*pow(d[j], 8.0) + 115.0/12.0*pow(d[j], 9.0);

             } else if( (d[j] >= 2.0) && (d[j] < 3.0) ){

               W44[j] = 5913.0 - 89235.0/4.0*d[j] + 297585.0/8.0*pow(d[j], 2.0) - 143895.0/4.0*pow( d[j], 3.0) + 177871.0/8.0*pow( d[j], 4.0) - 54641.0/6.0*pow(d[j], 5.0) + 19775.0/8.0*pow(d[j], 6.0) - 1715.0/4.0*pow(d[j], 7.0) + 345.0/8.0*pow(d[j], 8.0) - 23.0/12.0*pow(d[j], 9.0);

             } else{

               W44[j] = 0;
             }

            }


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



