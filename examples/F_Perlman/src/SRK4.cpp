#include "SRK4.H"
#include <iostream>
#include <vector>
#include "particle.hpp"
#include <cmath>
#include "interp.H"

State::State(const vector<particle> a,  const double a_dt,const double a_hp,  const double a_dx, const double domain, const double w)
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
     double z,x,y,omega;
     double c = 1.5;
     long double sumpsi = 0;
     double sumUg = 0, sumUg1 = 0;


     interpol.W44(vort, a_State.X, a_State.m_dx, a_State.hp, Np);
     double h2 = a_State.L/(n+1); 
     int px, py;
     Point i_closest;
     psi.setToZero();
     a_State.corrVort_error.clear();
     a_State.corrRelVort.clear();
     for(int i = 0; i < Np; i++){

       px = static_cast<int> ( round(a_State.X.at(i).x[0] / a_State.m_dx) );
       py = static_cast<int> ( round(a_State.X.at(i).x[1] / a_State.m_dx) );

       i_closest = Point({px, py} );
     
       omega = pow(1 - ( pow(a_State.X.at(i).x[1]-1.5, 2.0)+pow(a_State.X.at(i).x[0]-1.5, 2.0)), 7.0);
       a_State.corrVort_error.push_back( abs(*vort.data(i_closest)-omega) );
       a_State.corrRelVort.push_back( abs(*vort.data(i_closest)-omega) / abs(omega) ); 

     }

     vort.copyTo(psi, B, B); hockney.convolve(psi); 
    

     for(int i = 0; i < Np; i++){

       px = static_cast<int> ( round(a_State.X.at(i).x[0] / a_State.m_dx) );
       py = static_cast<int> ( round(a_State.X.at(i).x[1] / a_State.m_dx) );

       i_closest = Point({px, py} );


     
     }
     for( int dir = 0; dir < DIM; dir++){
        m_derivative[dir] = Stencil<double>::Derivative(1,dir,2);
     }
     
     BoxData<double> Ug(B); Ug.setToZero();
     Ug  = m_derivative[1](psi, 1.0/a_State.m_dx);
     BoxData<double> Vg(B); Vg.setToZero();
     Vg = m_derivative[0](psi, 1.0/a_State.m_dx);
     Vg *= -1;
     Point r;


     for( auto it = B2.begin(); it != B2.end(); it++){

       r = *it;
       z = pow( pow(r[0]*a_State.m_dx - 1.5, 2.0) + pow(r[1]*a_State.m_dx -1.5, 2.0), 0.5 );

       x = r[0]*a_State.m_dx;
       y = r[1]*a_State.m_dx;
       if( z<= 1){

         sumUg += (pow( *Ug.data(r) -  (y-c)/(16*pow(z, 2.0))*(1 - pow(1 - pow(z, 2.0), 8.0)), 2.0) + pow( *Vg.data(r) +  (x-c)/(16*pow(z, 2.0))*(1 - pow(1 - pow(z, 2.0), 8.0)), 2.0))*pow(a_State.m_dx, 2.0);
         sumUg1 += (abs( *Ug.data(r) -  (y-c)/(16*pow(z, 2.0))*(1 - pow(1 - pow(z, 2.0), 8.0))) + abs( *Vg.data(r) +  (x-c)/(16*pow(z, 2.0))*(1 - pow(1 - pow(z, 2.0), 8.0))))*pow(a_State.m_dx, 2.0);

       }else{

         sumUg += (pow( *Ug.data(r) - (y-c)/(16*pow(z, 2.0)), 2.0) + pow( *Vg.data(r) + (x-c)/(16*pow(z, 2.0)), 2.0))*pow(a_State.m_dx, 2.0);
         sumUg1 += (abs( *Ug.data(r) - (y-c)/(16*pow(z, 2.0))) + abs( *Vg.data(r) + (x-c)/(16*pow(z, 2.0))))*pow(a_State.m_dx, 2.0);


       }

     }

       
       sumUg = pow( sumUg, 0.5);
       a_State.sumUg = sumUg;
       a_State.sumUg1 = sumUg1;


     

     interpol.W44p( Ug, u, a_State.X, a_State.hp, a_State.m_dx, Np);
     interpol.W44p( Vg, v, a_State.X, a_State.hp, a_State.m_dx, Np);

    
}


void F::operator()( DX& a_DX, double a_time, double a_dt, State& a_State){
    
     particle p;
     vector<particle> Y, vel;
     vector<double> u,v;
     vector<particle> Pg;
     int n = a_State.L/a_State.m_dx; //number of grid points - 1
     Box Bg(Point::Zeros(), Point(n,n) );   
     Box B2( Point({1,1}), Point({n-1, n-1}) );
     BoxData<double> psi(Bg);
     BoxData<double> Vort_error(Bg);
     BoxData<double> Ug_error(B2);
     BoxData<double> Vg_error(B2);
     interp interpol;
     int  Np = a_State.X.size();
     int M = log2(n+1);
     double wg = 0;
     Hockney hockney(a_State.m_dx, M);
     double tol = pow(10.0, -10.0);
     double sumVort = 0;
     double sumU = 0;
     double sumUg = 0, sumUg1 = 0, sumU1 = 0, sumVort1 = 0;
     double sumpsi = 0;
     double z, z2;
     double c = 1.5;
     int step = .16875/a_dt;
     int t_pt;
     double x,y;
     t_pt = (a_time-a_dt)/a_dt;
     double omega;
     int px, py;
     Point i_closest;     
     
     double hp = a_State.hp;
     double h = a_State.m_dx;
     double h2 =a_State.L/(n-1);
     a_State.corrVort_error.clear();
     Vort_error.setToZero(); Ug_error.setToZero();
     Vg_error.setToZero();
     a_State.corrRelVort.clear();
     Point r;
     Stencil<double> m_derivative[DIM];
 
     for(unsigned int i = 0; i < a_DX.X2.size(); i++)
     {
          p.x[0] = a_DX.X2.at(i).x[0] + a_State.X.at(i).x[0];
          p.x[1] = a_DX.X2.at(i).x[1] + a_State.X.at(i).x[1];
          p.strength = a_State.X.at(i).strength;

          Y.push_back(p);

     }

     //Interpolate to grid 
     interpol.W44(psi, a_State.X, a_State.m_dx, a_State.hp, Np);


     for(int i = 0; i < Np; i++){

       px = static_cast<int> ( round(a_State.X.at(i).x[0] / a_State.m_dx) );
       py = static_cast<int> ( round(a_State.X.at(i).x[1] / a_State.m_dx) );

       i_closest = Point({px, py} );

       omega = pow(1 - ( pow(a_State.X.at(i).x[1]-1.5, 2.0)+pow(a_State.X.at(i).x[0]-1.5, 2.0)), 7.0);
       a_State.corrVort_error.push_back( abs(*psi.data(i_closest)-omega) );

       a_State.corrRelVort.push_back(abs(*psi.data(i_closest)-omega)/abs(omega) );

     }



     //Calculate voriticty summation over grid and error
     for( auto it = Bg.begin(); it != Bg.end(); it++){
       r = *it;

         

       wg += *psi.data(r)*pow(a_State.m_dx, 2.0);

       z = pow( pow(r[0]*a_State.m_dx - 1.5, 2.0) + pow(r[1]*a_State.m_dx -1.5, 2.0), 0.5 );
       
       if( z <= 1){
          Vort_error(r) += abs(*psi.data(r) - pow(1-pow(z, 2.0), 7.0));
          sumVort  +=  pow( *psi.data(r) - pow(1-pow(z, 2.0), 7.0), 2.0)*pow(a_State.m_dx, 2.0);

	  sumVort1 +=  abs( *psi.data(r) - pow(1-pow(z, 2.0), 7.0))*pow(a_State.m_dx, 2.0);

        } else{
   
          Vort_error(r) += abs(*psi.data(r));

          sumVort += pow( *psi.data(r), 2.0)*pow(a_State.m_dx, 2.0);
          sumVort1 += abs( *psi.data(r))*pow(a_State.m_dx, 2.0);

        }


     }


     //Check conservation
     if( abs(a_State.wp - wg) > pow(10.0, -9.0)){

       cout << "conservation of vorticity not guaranteed after interpolation to grid in RK4 step at time " << a_time << " with error " << abs(wg-a_State.wp) << endl;

     }

     

   if( t_pt%step == 0){
     //Use Hockney to calculate streamfunction
     #ifdef PR_HDF5
      HDF5Handler h4;
      h4.writePatch(a_State.m_dx,psi, "vortg%i.hdf5", t_pt/step);
      h4.writePatch(a_State.m_dx, Vort_error, "vort_error%i.hdf5", t_pt/step);
     #endif
    }
  
     hockney.convolve(psi);

     //Determine velocity on the grid

     for( int dir = 0; dir < DIM; dir++){
        m_derivative[dir] = Stencil<double>::Derivative(1,dir,2);
     }

     BoxData<double> Ug(Bg);
     Ug.setToZero();
     Ug  = m_derivative[1](psi, 1.0/a_State.m_dx);
     BoxData<double> Vg(Bg);
     Vg.setToZero();
     Vg = m_derivative[0](psi, 1.0/a_State.m_dx );
     Vg *= -1;

//       cout << a_time << " " << a_dt << " " << t_pt << endl;
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

     //Interpolate from grid to particles
     interpol.W44p( Ug, u, Y, a_State.hp, a_State.m_dx, Np);
     interpol.W44p( Vg, v, Y, a_State.hp, a_State.m_dx, Np);
     

     //velocity of particles error
     for(int i = 0; i < Np; i++){

       	     
       z =( pow( pow((a_State.X.at(i).x[0] - c),2.0) + pow((a_State.X.at(i).x[1] - c),2.0), 0.5 ));

       if( z <= 1){

         sumU += ( pow( u[i] + (a_State.X.at(i).x[1]-c )/(16*pow(z, 2.0) )*(1-pow(1-pow(z, 2.0), 8.0 )), 2.0 ) + pow( v[i]-(a_State.X.at(i).x[0]-c)/(16*pow(z, 2.0) )*(1-pow(1-pow(z, 2.0), 8.0 )), 2.0 ))*pow(hp, 2.0);

	 sumU1 +=( abs( u[i] + (a_State.X.at(i).x[1]-c )/(16*pow(z, 2.0) )*(1-pow(1-pow(z, 2.0), 8.0 ))) + abs( v[i]-(a_State.X.at(i).x[0]-c)/(16*pow(z, 2.0) )*(1-pow(1-pow(z, 2.0), 8.0 ))))*pow(hp, 2.0);

       }else{

         sumU += ( pow( u[i] + (a_State.X.at(i).x[1]-c)/(16*pow(z, 2.0) ), 2.0 ) + pow( v[i] -(a_State.X.at(i).x[0] - c)/(16*pow(z, 2.0) ), 2.0 ) )*pow(hp, 2.0);
         sumU1 += ( abs( u[i] + (a_State.X.at(i).x[1]-c)/(16*pow(z, 2.0) )) + abs( v[i] -(a_State.X.at(i).x[0] - c)/(16*pow(z, 2.0) ) ) )*pow(hp, 2.0);


    }

   }

     for( auto it = B2.begin(); it != B2.end(); it++){

       r = *it;
       z = pow( pow(r[0]*a_State.m_dx - 1.5, 2.0) + pow(r[1]*a_State.m_dx -1.5, 2.0), 0.5 );

       x = r[0]*a_State.m_dx;
       y = r[1]*a_State.m_dx;
       if( z<= 1){
       
	 Ug_error(r) += abs( *Ug.data(r) - (y-c)/(16.9*pow(z, 2.0))*(1 - pow(1 - pow(z, 2.0), 8.0)) );
	 Vg_error(r) += abs(*Vg.data(r) - (x-c)/(16.0*pow(z, 2.0))*(1 - pow(1 - pow(z, 2.0), 8.0)) ); 
         sumUg += (pow( *Ug.data(r) -  (y-c)/(16.0*pow(z, 2.0))*(1 - pow(1 - pow(z, 2.0), 8.0)), 2.0) + pow( *Vg.data(r) +  (x-c)/(16.0*pow(z, 2.0))*(1 - pow(1 - pow(z, 2.0), 8.0)), 2.0))*pow(a_State.m_dx, 2.0);
         sumUg1 += (abs( *Ug.data(r) -  (y-c)/(16.0*pow(z, 2.0))*(1 - pow(1 - pow(z, 2.0), 8.0))) + abs( *Vg.data(r) +  (x-c)/(16.0*pow(z, 2.0))*(1 - pow(1 - pow(z, 2.0), 8.0))))*pow(a_State.m_dx, 2.0);

       }else{

         Ug_error(r) += abs( *Ug.data(r) - (y-c)/(16.0*pow(z, 2.0)) );
         Vg_error(r) += abs( *Vg.data(r) + (x-c)/(16.0*pow(z, 2.0)) );
         sumUg += (pow( *Ug.data(r) - (y-c)/(16.0*pow(z, 2.0)), 2.0) + pow( *Vg.data(r) + (x-c)/(16.0*pow(z, 2.0)), 2.0))*pow(a_State.m_dx, 2.0);
         sumUg1 += (abs( *Ug.data(r) - (y-c)/(16.0*pow(z, 2.0))) + abs( *Vg.data(r) + (x-c)/(16.0*pow(z, 2.0))))*pow(a_State.m_dx, 2.0);

       }

  

     }


     if( t_pt%step == 0){
       #ifdef PR_HDF5
       HDF5Handler h6;
       h6.writePatch(a_State.m_dx,Ug_error, "Ugerror%i.hdf5", t_pt/step);
       h6.writePatch(a_State.m_dx,Vg_error, "Vgerror%i.hdf5", t_pt/step);
       #endif


    }



     sumVort = pow( sumVort, 0.5); 
     sumU = pow( sumU, 0.5);
     sumUg = pow(sumUg, 0.5);
     a_State.sumU = sumU;
     a_State.sumUg = sumUg;
     a_State.sumVort = sumVort;
     a_State.sumU1 = sumU1;
     a_State.sumUg1 = sumUg1;
     a_State.sumVort1 = sumVort1;


//     Vort_error.copyTo(a_State.Vort_error);

     //generate updated struct for k
     for( int i = 0;  i < Np; i++){

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
	      
        p.x[0] = r[0]*hp; p.x[1] = r[1]*hp;
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
