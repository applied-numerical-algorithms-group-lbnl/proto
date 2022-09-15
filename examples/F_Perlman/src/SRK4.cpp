#include "SRK4.H"
#include <iostream>
#include <vector>
#include "particle.hpp"
#include <cmath>
#include "interp.H"
#include "Proto_Timer.H"
#include <array>
#include "writers.h"

State::State(const vector<particle> a,  const double a_dt,const double a_hp,  const double a_dx, const double domain, const double w, const vector<double> lam)
{

        m_dt = a_dt;
	m_dx = a_dx;
        hp = a_hp;
	particle p;
	L = domain;
	wp = w;
	

        for(unsigned int i = 0; i < a.size(); i++){
         p.x[0] = a.at(i).x[0]; p.x[1] = a.at(i).x[1];
	 p.strength = a.at(i).strength;

	 lambda.push_back(lam.at(i));
	 X.push_back(p);
	
	}

;

    
};
void 
State::increment(const DX& a_DX)
{  
 
 for(unsigned int i = 0; i < X.size(); i++)
  { 
    X.at(i).x[0] += a_DX.X2.at(i).x[0];
    X.at(i).x[1] += a_DX.X2.at(i).x[1];
     
  }

  
};
void 
DX::increment(double& a_weight,const DX& a_incr)
{
   PR_TIMERS("main::RK4::increment");
   for(unsigned int i = 0; i < X2.size(); i++){
     X2.at(i).x[0] += a_weight*a_incr.X2.at(i).x[0];
     X2.at(i).x[1] += a_weight*a_incr.X2.at(i).x[1];
   }


};

void 
DX::init(State& a_State)
{
    particle p;
    
    X2 = a_State.X;

    for(unsigned int i = 0; i < a_State.X.size(); i++){

            X2.at(i).x[0] = 0; 
	    X2.at(i).x[1] = 0;
  
    }

   
};

void 
DX::operator*=(const double& a_weight)
{

    for(unsigned int i = 0; i < X2.size(); i++){
	X2.at(i).x[0] *= a_weight;
        X2.at(i).x[1] *= a_weight;

    }


};

void State::getVelocity(BoxData<double>& vort,vector<double>& errorVg, vector<double>& errorpsi, vector<double>& u, vector<double>& v, const int Np, int step, State& a_State ){

     int n = a_State.L/a_State.m_dx; //number of grid points - 1
     Box B(Point::Zeros(), Point(n,n) );
     Box B2( Point({1,1}), Point({n-1, n-1}) );
     int np = a_State.L/a_State.hp;
     BoxData<double> psi(B);
     BoxData<double> Vf(B);
     BoxData<double> Uf(B);
     BoxData<double> psi_exact(B);
     BoxData<double> psi_error(B);
     interp interpol;
     int M = log2(n+1);
     Stencil<double> m_derivative[DIM];
     Hockney hockney(a_State.m_dx, M);
     double z,x,y,omega, ux, uy, z2;
     double c = 2.0;
     long double sumpsi = 0;
     double sumUg = 0, sumUg1 = 0;

     //Interpolate from grid to particles
     interpol.W44(vort, a_State.X, a_State.m_dx, a_State.hp, Np);

     double h2 = pow(a_State.m_dx, 2.0); 
     int px, py;

     Point i_closest;
     psi.setToZero();
     a_State.radi.clear();
     a_State.corrVort_error.clear();
     a_State.corrRelVort.clear();

     //Calculate relative and exact error in vorticity
     for(int i = 0; i < Np; i++){

       px = static_cast<int> ( round((a_State.X.at(i).x[0]+c) / a_State.m_dx) );
       py = static_cast<int> ( round((a_State.X.at(i).x[1]+c) / a_State.m_dx) );

       i_closest = Point({px, py} );

     
       omega = pow(1 - ( pow(a_State.X.at(i).x[1], 2.0)+pow(a_State.X.at(i).x[0], 2.0)), 7.0);
       a_State.corrVort_error.push_back( abs(*vort.data(i_closest)-omega) );
       a_State.corrRelVort.push_back( abs(*vort.data(i_closest)-omega) / abs(omega) ); 
       a_State.radi.push_back(pow(pow(a_State.X.at(i).x[1], 2.0)+pow(a_State.X.at(i).x[0], 2.0), 0.5));

     }

     vort.copyTo(psi, B, B); hockney.convolve(psi); 
    
     //Calculate velocity from Hockney results
      for( int dir = 0; dir < DIM; dir++){
        m_derivative[dir] = Stencil<double>::Derivative(1,dir,2);
     }
     
     BoxData<double> Ug(B); Ug.setToZero();
     Ug  = m_derivative[1](psi, 1.0/a_State.m_dx);
     BoxData<double> Vg(B); Vg.setToZero();
     Vg = m_derivative[0](psi, 1.0/a_State.m_dx);
     Vg *= -1;
     Point r;


     //Error in velocity on the grid
     for( auto it = B2.begin(); it != B2.end(); it++){

       r = *it;
       z2 = ( pow(r[0]*a_State.m_dx - c, 2.0) + pow(r[1]*a_State.m_dx -c, 2.0) );

       x = r[0]*a_State.m_dx-c;
       y = r[1]*a_State.m_dx-c;

       ux = y/(16.0*z2)*(1.0-pow(1.0-z2, 8.0));
       uy = -x/(16.0*z2)*(1.0-pow(1.0-z2, 8.0));
       if( z<= 1){

         sumUg += (pow( *Ug.data(r) -  ux, 2.0) + pow( *Vg.data(r) - uy, 2.0))*h2;
         sumUg1 += (abs( *Ug.data(r) -  ux) + abs( *Vg.data(r) - uy))*h2;

       }else{

         sumUg += (pow( *Ug.data(r) - (y)/(16*z2), 2.0) + pow( *Vg.data(r) + (x)/(16*z2), 2.0))*h2;
         sumUg1 += (abs( *Ug.data(r) - (y)/(16*z2)) + abs( *Vg.data(r) + (x)/(16*z2)))*h2;


       }

     }

       
       sumUg = pow( sumUg, 0.5);
       a_State.sumUg = sumUg;
       a_State.sumUg1 = sumUg1;


     //Interpolate to get particle velocity
     interpol.W44p( Ug, u, a_State.X, a_State.hp, a_State.m_dx, Np);
     interpol.W44p( Vg, v, a_State.X, a_State.hp, a_State.m_dx, Np);
    
}


void F::operator()( DX& a_DX, double a_time, double a_dt, State& a_State){
    
     particle p;
     particle lam;
     vector<particle> Y, vel;
     vector<double> u,v;
     vector<particle> Pg;
     vector<array<double, DIM>> corr;
     vector<particle> lambda_p;

     int n = a_State.L/a_State.m_dx; //number of grid points - 1
     int Np = a_State.X.size();
     int M = log2(n+1);
     int step = .16875/a_dt, t_pt;
     int px, py;

     Box Bg(Point::Zeros(), Point(n,n) );   
     Box B2( Point({2,2}), Point({n-2, n-2}) );
     BoxData<double> psi(Bg);
     BoxData<double> lambda_g(Bg);
     BoxData<double> Vort_error(Bg);
     BoxData<double> Ug_error(B2);
     BoxData<double> Vg_error(B2);
     BoxData<double> lamb_g(Bg);
     interp interpol;
     Hockney hockney(a_State.m_dx, M);

     double tol = pow(10.0, -10.0);
     double sumVort = 0, sumVortm,  sumU = 0;
     double sumUg = 0, sumUg1 = 0, sumU1 = 0, sumVort1 = 0;
     double sumpsi = 0, wg = 0;
     double z, z2;
     double c = 2.0;
     double x,y;
     double omega, ux, uy, w_o;
     double hp = a_State.hp;
     double h = a_State.m_dx;
     double h2 = pow(a_State.m_dx, 2.0);
     Point i_closest; Point r;

     const char *const varnames[] = {"LambdaVort_Error", "Vorticity_Error"};     
     string filename = "vorticity_lambda_corr";
  
     t_pt = (a_time-a_dt)/a_dt;
     a_State.corrVort_error.clear();
     Vort_error.setToZero(); Ug_error.setToZero();
     Vg_error.setToZero(); lambda_g.setToZero();
     a_State.corrRelVort.clear();
     Stencil<double> m_derivative[DIM];
 
     //Update particle and lambda position with RK4 step 
     for(unsigned int i = 0; i < a_DX.X2.size(); i++){

          p.x[0] = a_DX.X2.at(i).x[0] + a_State.X.at(i).x[0];
          p.x[1] = a_DX.X2.at(i).x[1] + a_State.X.at(i).x[1];
          p.strength = a_State.X.at(i).strength;

          lam.x[0] = a_State.X.at(i).x[0];
          lam.x[1] = a_State.X.at(i).x[1];
          lam.strength = 1;

	  lambda_p.push_back(lam);
          Y.push_back(p);

     }


     //Interpolate to grid
     { PR_TIMERS("main::RK4::interp"); 
     interpol.W44(psi, a_State.X, a_State.m_dx, a_State.hp, Np);
     interpol.W44(lambda_g, lambda_p, a_State.m_dx, a_State.hp, Np);
     }


     //Calculate error in vorticity field
     for(int i = 0; i < Np; i++){
       px = static_cast<int> ( round((a_State.X.at(i).x[0]+c) / a_State.m_dx) );
       py = static_cast<int> ( round((a_State.X.at(i).x[1]+c) / a_State.m_dx) );

       i_closest = Point({px, py} );

       omega = pow(1 - ( pow(a_State.X.at(i).x[1], 2.0)+pow(a_State.X.at(i).x[0], 2.0)), 7.0);
       a_State.corrVort_error.push_back( abs(*psi.data(i_closest)-omega) );

       a_State.corrRelVort.push_back(abs(*psi.data(i_closest)-omega)/abs(omega) );

     }



     array<double,DIM> temp;
     BoxData<double> lambda_error(Bg); lambda_error.setToZero();
     //Calculate voriticty summation over grid and error
     for( auto it = Bg.begin(); it != Bg.end(); it++){
       r = *it;


       wg += *psi.data(r)*pow(a_State.m_dx, 2.0);

       z = pow( pow(r[0]*a_State.m_dx - c, 2.0) + pow(r[1]*a_State.m_dx -c, 2.0), 0.5 );
      
       omega = pow(1-pow(z, 2.0), 7.0); 

       if( z <= 1){

	  if(abs(*psi.data(r) - omega) > sumVortm){

            sumVortm = abs(*psi.data(r) - omega);

	  }
          Vort_error(r) += abs(*psi.data(r) - omega);
          sumVort  +=  pow( *psi.data(r) - omega, 2.0)*h2;

	  sumVort1 +=  abs( *psi.data(r) - omega)*h2;
	  temp[0] =abs(*psi.data(r) - omega);

        } else{
 
	  if(abs(*psi.data(r) ) > sumVortm){

            sumVortm = abs(*psi.data(r) );

          }

	  temp[0] = abs(*psi.data(r) );
	
          Vort_error(r) += abs(*psi.data(r));

          sumVort += pow( *psi.data(r), 2.0)*h2;
          sumVort1 += abs( *psi.data(r))*h2;

        }

       temp[1] = abs(abs(*lambda_g.data(r) -1)*(*psi.data(r)));

        lamb_g(r) = (*psi.data(r))*(*lambda_g.data(r));
        lambda_error(r) += abs(abs(*lambda_g.data(r) -1)*(*psi.data(r)));
        corr.push_back(temp);

     }

     PWrite<DIM>(corr, varnames, filename, t_pt);

    //Check conservation
     if( abs(a_State.wp - wg) > pow(10.0, -9.0)){

       cout << "conservation of vorticity not guaranteed after interpolation to grid in RK4 step at time " << a_time << " with error " << abs(wg-a_State.wp) << endl;

     }

     

   if( t_pt%step == 0){

     //Use Hockney to calculate streamfunction
     #ifdef PR_HDF5
      HDF5Handler h4;
      h4.writePatch(a_State.m_dx, lambda_g, "Lambdag%i.hdf5", t_pt/step);
      h4.writePatch(a_State.m_dx,psi, "vortg%i.hdf5", t_pt/step);
      h4.writePatch(a_State.m_dx, Vort_error, "vort_error%i.hdf5", t_pt/step);
      h4.writePatch(a_State.m_dx, lambda_error, "LambdaVort_error%i.hdf5", t_pt/step);
     #endif
    }
  
     hockney.convolve(psi);

     //Determine velocity on the grid
     for( int dir = 0; dir < DIM; dir++){
        m_derivative[dir] = Stencil<double>::Derivative(1,dir,4);
     }
     

     BoxData<double> Ug(Bg);
     BoxData<double> Vg(Bg);    
     {
     Ug.setToZero();
     Ug  = m_derivative[1](psi, 1.0/a_State.m_dx);
     Vg.setToZero();
     Vg = m_derivative[0](psi, 1.0/a_State.m_dx );
     Vg *= -1;
     }
     
    //Exact Velocity on th grid
    //for( auto it = Bg.begin(); it != Bg.end(); it++){
    //r = it.operator*();
    //z = ( pow( pow((r[0]*a_State.m_dx - c),2.0) + pow((r[1]*a_State.m_dx - c),2.0), 0.5 ));
    //if( z <= 1){
    //Ug(r) += 1.0/(16.0*pow(z, 2.0) )*(1 - pow( (1 - pow(z, 2.0) ), 8.0) )*(r[1]*a_State.m_dx - c);
    //Vg(r) += -1.0/(16.0*pow(z, 2.0) )*(1 - pow( (1 - pow(z, 2.0) ), 8.0) )*(r[0]*a_State.m_dx - c);
    //} else{
    //z2 = pow((r[0]*a_State.m_dx - c),2.0) + pow((r[1]*a_State.m_dx - c),2.0);
    //Ug(r) += 1.0/(16.0*z2) *(r[1]*a_State.m_dx - c);
    //Vg(r) += -1.0/(16.0*z2) *(r[0]*a_State.m_dx - c);
    //}
    //}


    if( t_pt%step == 0){
       #ifdef PR_HDF5
       HDF5Handler h5;
       h5.writePatch(a_State.m_dx,Ug, "Ug%i.hdf5", t_pt/step);
       h5.writePatch(a_State.m_dx,Vg, "Vg%i.hdf5", t_pt/step);
       #endif

       #ifdef PR_HDF5
       HDF5Handler h3;
       h3.writePatch(a_State.m_dx,psi, "psi%i.hdf5", t_pt/step);
       #endif

    }
 


     { //Interpolate from grid to particles

     PR_TIMERS("main::RK4::interp_grid");
     interpol.W44p( Ug, u, Y, a_State.hp, a_State.m_dx, Np);
     interpol.W44p( Vg, v, Y, a_State.hp, a_State.m_dx, Np);

     }


     //Calculate error on grid
     for( auto it = B2.begin(); it != B2.end(); it++){

      PR_TIMERS("main::RK4::UG_error"); 
      r = *it;
      z2 = ( pow(r[0]*a_State.m_dx -c, 2.0) + pow(r[1]*a_State.m_dx -c, 2.0) );
     

      x = r[0]*a_State.m_dx-c;
      y = r[1]*a_State.m_dx-c;

      ux = (y)/(16.0*z2)*(1 - pow(1 - z2, 8.0));
      uy = (-1.0*x)/(16.0*z2)*(1 - pow(1 - z2, 8.0));

      if( z2 <= 1){
       
         Ug_error(r) += abs( *Ug.data(r) - ux );
         Vg_error(r) += abs(*Vg.data(r) - uy ); 
         sumUg += (pow( *Ug.data(r) -  ux, 2.0) + pow( *Vg.data(r) - uy, 2.0))*h2;
         sumUg1 += (abs( *Ug.data(r) -  ux) + abs( *Vg.data(r) - uy))*h2;

       }else{

         Ug_error(r) += abs( *Ug.data(r) - (y)/(16.0*z2) );
         Vg_error(r) += abs( *Vg.data(r) + (x)/(16.0*z2) );
         sumUg += (pow( *Ug.data(r) - (y)/(16.0*z2), 2.0) + pow( *Vg.data(r) + (x)/(16.0*z2), 2.0))*h2;
         sumUg1 += (abs( *Ug.data(r) - (y)/(16.0*z2) + abs( *Vg.data(r) + (x)/(16.0*z2))))*h2;

       }

     }

//     if( t_pt%step == 0){
//       #ifdef PR_HDF5
//       HDF5Handler h6;
//       h6.writePatch(a_State.m_dx,Ug_error, "Ugerror%i.hdf5", t_pt/step);
//       h6.writePatch(a_State.m_dx,Vg_error, "Vgerror%i.hdf5", t_pt/step);
//       #endif
//    }



     sumVort = pow( sumVort, 0.5); 
     sumU = pow( sumU, 0.5);
     sumUg = pow(sumUg, 0.5);
     a_State.sumU = sumU;
     a_State.sumUg = sumUg;
     a_State.sumVort = sumVort;
     a_State.sumVortm = sumVortm;
     a_State.sumU1 = sumU1;
     a_State.sumUg1 = sumUg1;
     a_State.sumVort1 = sumVort1;



     //generate updated info for k
     for( int i = 0;  i < Np; i++){

       //Stuff for exact particle velocity for debugging	     
       
       //z =  pow( pow(Y.at(i).x[1] ,2.0) + pow(Y.at(i).x[0] ,2.0), 0.5 );
       //if( (abs(Y.at(i).x[1]) < pow(10.0, -12.0)) && (abs(Y.at(i).x[0]) < pow(10.0, -12.0)) ){
       //p.x[0] = 0.0;
       //p.x[1] = 0.0;
       //}else if( z <= 1){
       //p.x[0] = (1.0/(16.0*pow(z, 2.0) )*(1 - pow( (1 - pow(z, 2.0) ), 8.0) )*(Y.at(i).x[1] ))*a_dt;
       //p.x[1]  = (-1.0/(16.0*pow(z, 2.0) )*(1 - pow( (1 - pow(z, 2.0) ), 8.0) )*(Y.at(i).x[0] ))*a_dt;
       //} else{
       //z2 = pow((Y.at(i).x[1] ),2.0) + pow((Y.at(i).x[0] ),2.0);
       //p.x[0] = ( 1.0/(16.0*z2) *(Y.at(i).x[1] ))*a_dt;
       //p.x[1]  = ( -1.0/(16.0*z2) *(Y.at(i).x[0] ))*a_dt;
       //}
     
       p.x[0] = a_dt*u[i]; p.x[1] =a_dt*v[i];
       p.strength = 0;

	vel.push_back(p);
     
     }


     a_DX.X2 = vel;
    
  
}

void State::remap(){


    interp interpol;
    int n = L/hp;
    int Np = X.size();
    Box Bp(Point::Zeros(), Point(n,n) );
    BoxData<double> strength_remap(Bp);
    Point r;
    double strength;
    double c = 2.0;
    double wc = 0;
    double k = 0;
    double wp2 = 0;


    //Deposit on grid
    interpol.W44(strength_remap,X, hp, hp, Np);

    X.clear();
    particle p;

    // Generate new particle vector
    for(auto it = Bp.begin(); it != Bp.end(); it++){
      
      r = *it;

      strength = *strength_remap.data(r);
  
      if( abs(strength) > pow(10, -9.0)){
	      
        p.x[0] = r[0]*hp-c; p.x[1] = r[1]*hp-c;
	p.strength = strength;
        k++;     

	X.push_back(p);
      }

    }


    //Check conservation
    for( auto it = Bp.begin(); it != Bp.end(); it++){
       r = *it;

       wc += *strength_remap.data(r)*pow(hp, 2.0);
    }
    

    if( abs(wp - wc) > pow(10.0, -14.0)){

       cout << "conservation of vorticity not guaranteed in remap " << abs(wc-wp) << endl;

    }

 

   
}
