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
#include "Proto_WriteBoxData.H"
#include "Proto_Timer.H"
#include <array>
#include "Proto_Timer.H"
#define NUMCOMPS DIM+2

using namespace std;
using namespace Proto;



/***/
int main(int argc, char* argv[])
{

    double error32 = 0, error64 = 0, error128 = 0, error64_2 = 0;
    int step = 88;
    int n = 31, n2 = 63, n3 = 127; //number of grid points - 1
    Box B32( Point({2,2}), Point({n-2, n-2}) );
    Box B32_2( Point{3,3}, Point{n-3, n-3});
    Box B64( Point({2,2}), Point({n2-2, n2-2}) );
    Box B128( Point({2,2}), Point({n3-2, n3-2}) );
    Box Bsten(Point({-1,-1}), Point({1,1}) );
    Box Bsten_2(Point({-2,-2}), Point({2,2}) );
    BoxData<double> Ug32(B32);
    BoxData<double> Vg32(B32);
    BoxData<double> Ug64(B64);
    BoxData<double> Vg64(B64);
    BoxData<double> Ug128(B128);
    BoxData<double> Vg128(B128);
    BoxData<double> Ug128to64(B64);
    BoxData<double> Vg128to64(B64);
    BoxData<double> Ug128to64_r(B64);
    BoxData<double> Vg128to64_r(B64);
    BoxData<double> Ug64to32(B32);
    BoxData<double> Vg64to32(B32);
    BoxData<double> Ug128to32(B32);
    BoxData<double> Vg128to32(B32);
    BoxData<double> Ug64to32_r(B32);
    BoxData<double> Vg64to32_r(B32);
    BoxData<double> Ugrichard(B32);
    BoxData<double> Vgrichard(B32);
    BoxData<double> Ugint(B32);
    BoxData<double> Vgint(B32);
    vector<double> verror64, verror128, verror32, verror64_2, time;
    double dt = 0.16875;

    Point r; 
    Point p,k, k_to;
    Point i_to;
    double tmp;
    double h1 =1.5/n;
    double h2 = 1.5/n2;
    double h3 = 1.5/n3;
    double x_to[2];
    Ug32.setToZero();
    Ug64.setToZero();
    Ug128.setToZero();
    Vg32.setToZero();
    Vg64.setToZero();
    Vg128.setToZero();
    Ug128to64.setToZero();
    Vg128to64.setToZero();
    Ug128to64_r.setToZero();
    Vg128to64_r.setToZero();
    Ug64to32.setToZero();
    Vg64to32.setToZero();
    Ug128to32.setToZero();
    Vg128to32.setToZero();
    Ug64to32_r.setToZero();
    Vg64to32_r.setToZero();
    Ugint.setToZero();
    Vgint.setToZero();
    Ugrichard.setToZero();
    Vgrichard.setToZero();
    string path;
    error32 = 0;
    error64 = 0;
    error64_2 = 0;
    error128 = 0;
    for(int i = 0; i < step; i++){


	    path = "32/Ug_" + to_string(i) + ".txt";
	    ifstream Ug_32(path);
	    path = "64/Ug_" + to_string(i) + ".txt";
            ifstream Ug_64(path);
	    path = "128/Ug_" + to_string(i) + ".txt";
            ifstream Ug_128(path);

	    path = "32/Vg_" + to_string(i) + ".txt";
            ifstream Vg_32(path);
            path = "64/Vg_" + to_string(i) + ".txt";
            ifstream Vg_64(path);
            path = "128/Vg_" + to_string(i) + ".txt";
            ifstream Vg_128(path);



	    for( auto it = B32.begin(); it != B32.end(); it++){

		    r = *it;

                    Ug_32 >> tmp;
		    Ug32(r) += tmp;
                    Vg_32 >> tmp;
		    Vg32(r) += tmp;

	    }


           for( auto it = B64.begin(); it != B64.end(); it++){
                    r = *it;

                    Ug_64 >> tmp;
                    Ug64(r) += tmp;
                    Vg_64 >> tmp;
                    Vg64(r) += tmp;

            }

           for( auto it = B128.begin(); it != B128.end(); it++){
                    r = *it;

                    Ug_128 >> tmp;
                    Ug128(r) += tmp;
                    Vg_128 >> tmp;
                    Vg128(r) += tmp;

            }


            for( auto it = B64.begin(); it != B64.end(); it++){
                    r = *it;

		    if( B64.onBoundary(r) ){

		    }else{

			   p = Point({round(r[0]*h2/h3),round(r[1]*h2/h3)});

			   for(auto j = Bsten.begin(); j != Bsten.end(); j++){

			      i_to = p + *j;

			      Ug128to64(r) += Ug128(i_to)/9;
                              Vg128to64(r) += Vg128(i_to)/9;


			   }


		    

                         Ug128to64_r(r) += (16*Ug128to64(r) - Ug64(r)) /15;
		         Vg128to64_r(r) += (16*Vg128to64(r) - Vg64(r)) /15;
                      //   error64_2 = error64_2 + (pow(h2*(Ug128to64_r(r)-Ug64(r)), 2.0) + pow(h2*(Vg128to64_r(r)-Vg64(r)), 2.0) );
                      //   error128 = error128 + (pow(h2*(Ug128to64_r(r)-Ug128to64(r)), 2.0) + pow(h2*(Vg128to64_r(r)-Vg128to64(r)), 2.0) );

                    }



            }


	     for( auto it = B32_2.begin(); it != B32_2.end(); it++){
                    r = *it;

                    if( B32_2.onBoundary(r) ){

                    }else{

                           k = Point({round(r[0]*h1/h3),round(r[1]*h1/h3)});
                           for(auto j = Bsten_2.begin(); j != Bsten_2.end(); j++){

                              k_to = k + *j;
                              Ug128to32(r) += Ug128(k_to)/16;
                              Vg128to32(r) += Vg128(k_to)/16;



                           }

                    }

	     }

            for( auto it = B32.begin(); it != B32.end(); it++){
                    r = *it;

                    if( B32.onBoundary(r) ){

                    }else{

                           p = Point({round(r[0]*h1/h2),round(r[1]*h1/h2)});
                           k = Point({round(r[0]*h1/h3),round(r[1]*h1/h3)});
                           for(auto j = Bsten.begin(); j != Bsten.end(); j++){

                              i_to = p + *j;
                              k_to = k + *j;
                              Ug64to32(r) += Ug64(i_to)/9;
                              Vg64to32(r) += Vg64(i_to)/9;
			      Ugint(r)    += Ug128to64_r(i_to)/9;
			      Vgint(r)    += Vg128to64_r(i_to)/9;



                           }


                           Ug64to32_r(r) += (16*Ug64to32(r) - Ug32(r)) /15;
                           Vg64to32_r(r) += (16*Vg64to32(r) - Vg32(r)) /15;
			   Ugrichard(r)  += (32*Ugint(r) -    Ug64to32_r(r)) /31;
			   Vgrichard(r)  += (32*Vgint(r) -    Vg64to32_r(r)) /31;

                           error32 = error32 + (pow(h1*(Ugrichard(r)-Ug32(r)), 2.0) + pow(h1*(Vgrichard(r)-Vg32(r)), 2.0) );
		           error64 = error64 + (pow(h1*(Ugrichard(r)-Ug64to32(r)), 2.0) + pow(h1*(Vgrichard(r)-Vg64to32(r)), 2.0) );
			   error128 = error128 + (pow(h1*(Ugrichard(r)-Ug128to32(r)), 2.0) + pow(h1*(Vgrichard(r)-Vg128to32(r)), 2.0) );

            
	    
		    }  

		    
            }  

	    error32 = pow(error32, 0.5);
	    error64 = pow(error64, 0.5);
	  //  error64_2 = pow(error64_2, 0.5);
	    error128  = pow(error128, 0.5);

	    verror32.push_back(error32);
	    verror64.push_back(error64);
	 //   verror64_2.push_back(error64_2);
	    verror128.push_back(error128);
            time.push_back(dt*i);




           #ifdef PR_HDF5
           HDF5Handler h5;
           h5.writePatch(h2,Ug128to64, "Ug128to64_%i.hdf5", i);
           h5.writePatch(h2,Vg128to64, "Vg128to64_%i.hdf5", i);
	   h5.writePatch(h2,Ug128to64_r, "Ug128to64_r_%i.hdf5", i);
           h5.writePatch(h2,Vg128to64_r, "Vg128to64_r_%i.hdf5", i);
           h5.writePatch(h2,Ug64to32, "Ug64to32_%i.hdf5", i);
           h5.writePatch(h2,Vg64to32, "Vg64to32_%i.hdf5", i);
           h5.writePatch(h2,Ug64to32_r, "Ug64to32_r_%i.hdf5", i);
           h5.writePatch(h2,Vg64to32_r, "Vg64to32_r_%i.hdf5", i);
           h5.writePatch(h2,Ug128to64, "Ug128to64_%i.hdf5", i);
           h5.writePatch(h2,Vg128to64, "Vg128to64_%i.hdf5", i);
           h5.writePatch(h2,Ug128to32, "Ug128to32_%i.hdf5", i);
           h5.writePatch(h2,Vg128to32, "Vg128to32_%i.hdf5", i);
	   h5.writePatch(h2,Ugrichard, "Ugrichard_%i.hdf5", i);
           h5.writePatch(h2,Vgrichard, "Vgrichard_%i.hdf5", i);

           #endif


           Ug32.setToZero();
           Ug64.setToZero();
           Ug128.setToZero();
           Vg32.setToZero();
           Vg64.setToZero();
           Vg128.setToZero();
	   Ug128to64_r.setToZero();
           Vg128to64_r.setToZero();
           Ug128to64.setToZero();
           Vg128to64.setToZero();
           Ug64to32.setToZero();
           Vg64to32.setToZero();
           Ug64to32_r.setToZero();
           Vg64to32_r.setToZero();
           Ugrichard.setToZero();
           Vgrichard.setToZero();
           Ug128to32.setToZero();
           Vg128to32.setToZero();
           Ugint.setToZero();
           Vgint.setToZero();

           error32 = 0;
	   error64 = 0;
	   error128 = 0;
//	   error64_2 = 0;



    }

  string q,e;
  //Write curve files
  string errorfile2 = "errorVelocity32.curve";
  ofstream f1(errorfile2);
  string errorfile4 = "errorVelocity64.curve";
//  string errorfile5 = "errorVelocity64_2.curve";
  string errorfile7 = "errorVelocity128.curve";
  ofstream f2(errorfile4);
//  ofstream f3(errorfile5);
  ofstream f4(errorfile7);

  q = "# TIME";
  e = "# Velocity32";


  if( f1.is_open() ) {

      //format in scientific notation
       f1.setf(std::ios::scientific, std::ios::floatfield);
       f1.precision(4);

       f1 << q << endl;


       f1 << e << endl;


      for(unsigned int i = 0; i < time.size(); i++){

          f1 << time.at(i) << " " << verror32.at(i)  <<  endl;
      }

   }

   e = "# Velocity64";


  if( f2.is_open() ) {

      //format in scientific notation
       f2.setf(std::ios::scientific, std::ios::floatfield);
       f2.precision(4);

       f2 << q << endl;


       f2 << e << endl;


      for(unsigned int i = 0; i < time.size(); i++){

          f2 << time.at(i) << " " << verror64.at(i)  <<  endl;
      }

   }

//   e = "# Velocity64_2";


//  if( f3.is_open() ) {

      //format in scientific notation
//       f3.setf(std::ios::scientific, std::ios::floatfield);
//       f3.precision(4);

//       f3 << q << endl;


//       f3 << e << endl;


//      for(unsigned int i = 0; i < time.size(); i++){

//          f3 << time.at(i) << " " << verror64_2.at(i)  <<  endl;
//      }

//   }

    e = "# Velocity128";


  if( f4.is_open() ) {

      //format in scientific notation
       f4.setf(std::ios::scientific, std::ios::floatfield);
       f4.precision(4);

       f4 << q << endl;


       f4 << e << endl;


      for(unsigned int i = 0; i < time.size(); i++){

          f4 << time.at(i) << " " << verror128.at(i)  <<  endl;
      }

   }


    


}
