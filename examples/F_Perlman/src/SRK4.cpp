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
     long double z;
     long double sumpsi = 0;
     double sumUg = 0;
     interpol.W22(vort, a_State.X, a_State.m_dx, a_State.hp, Np);
    double h2 = a_State.L/(n+1); 
     psi.setToZero();

     vort.copyTo(psi, B, B); hockney.convolve(psi); 
    
   
     for( int dir = 0; dir < DIM; dir++){
        m_derivative[dir] = Stencil<double>::Derivative(1,dir,2);
     }
     
     BoxData<double> Ug(B); Ug.setToZero();
     Ug  = m_derivative[1](psi, 1.0/a_State.m_dx);
     BoxData<double> Vg(B); Vg.setToZero();
     Vg = m_derivative[0](psi, 1.0/a_State.m_dx);
     Vg *= -1;
     Point r;
     
//     Vg.copyTo(Vf, B, B); Ug.copyTo(Vf, B, B);



     psi_exact.setToZero(); psi_error.setToZero();

     for( auto it = B.begin(); it != B.end(); it++){
        r = it.operator*();


	z = pow( pow( (r[0]*a_State.m_dx - 1.5), 2.0) + pow( (r[1]*a_State.m_dx - 1.5), 2.0), 0.5 );
       
      	if( z <= 1){

           // sumUg += (pow(*Ug.data(r) - 0.5*(r[1]*a_State.m_dx - 1.5), 2.0)+pow(*Vg.data(r) + 0.5*(r[0]*a_State.m_dx-1.5), 2.0) )*pow(a_State.hp, 2.0);
	  
	    sumpsi += abs( (*psi.data(r) - 0.25*( pow(z, 2.0)) )+0.25)*pow(a_State.m_dx, 2.0);

	    psi_exact(r) +=  0.25*( pow(z, 2.0))-0.25;

	    

	    psi_error(r) = abs( *psi.data(r) - *psi_exact.data(r) );



        } else{

	  // sumUg += (pow(*Ug.data(r) - 0.5/pow(z, 2.0)*(r[1]*a_State.m_dx - 1.5), 2.0)+pow(*Vg.data(r) + 0.5/pow(z, 2.0)*(r[0]*a_State.m_dx-1.5), 2.0) )*pow(a_State.hp, 2.0);
           sumpsi += abs( (*psi.data(r) - 0.25*log( pow(z, 2.0) ) ))*pow( a_State.m_dx, 2.0);
   
	   psi_exact(r) += 0.25*log( pow(z, 2.0) );
	    
           psi_error(r) = abs( *psi.data(r) - *psi_exact.data(r) );
	}

//cout << r << " " << psi_exact(r) << endl;


       }
	 sumpsi = pow( sumpsi, 0.5);
   
   //  cout << sumpsi << endl;
     errorpsi.push_back(sumpsi);


      for( auto it = B2.begin(); it != B2.end(); it++){
        r = it.operator*();


        z = pow( pow( (r[0]*a_State.m_dx - 1.5), 2.0) + pow( (r[1]*a_State.m_dx - 1.5), 2.0), 0.5 );

        if( z <= 1){

            sumUg += (abs(*Ug.data(r) - 0.5*(r[1]*a_State.m_dx - 1.5))+abs(*Vg.data(r) + 0.5*(r[0]*a_State.m_dx-1.5)) )*pow(a_State.m_dx, 2.0);

        } else{

           sumUg += (abs(*Ug.data(r) - 0.5/pow(z, 2.0)*(r[1]*a_State.m_dx - 1.5))+abs(*Vg.data(r) + 0.5/pow(z, 2.0)*(r[0]*a_State.m_dx-1.5)) )*pow(a_State.m_dx, 2.0);
          
        }


       }
   
  
      errorpsi.push_back(sumpsi);
      errorVg.push_back(sumUg);

   //out << "interp" << endl;
     interpol.W22p( Ug, u, a_State.X, a_State.hp, a_State.m_dx, Np);
     interpol.W22p( Vg, v, a_State.X, a_State.hp, a_State.m_dx, Np);


//     #ifdef PR_HDF5
//     HDF5Handler h5;
//     h5.writePatch(a_State.m_dx,psi_exact, "psi_exact%i.hdf5", 0);
//     h5.writePatch(a_State.m_dx,psi_error, "psi_error%i.hdf5", step);
 //    #endif


    
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
     interp interpol;
     int  Np = a_State.X.size();
     int M = log2(n+1);
     double wg = 0;
     Hockney hockney(a_State.m_dx, M);
     double tol = pow(10.0, -10.0);
     double sumVort = 0;
     double sumU = 0;
     double sumUg = 0;
     double sumpsi = 0;
     double z, z2;
     int step = .16875/a_dt;
     int t_pt = (a_time-a_dt)/a_dt;

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
     interpol.W22(psi, a_State.X, a_State.m_dx, a_State.hp, Np);


     //Calculate voriticty summation over grid and error
     for( auto it = Bg.begin(); it != Bg.end(); it++){
       r = it.operator*();

       wg += *psi.data(r)*pow(a_State.m_dx, 2.0);

       if( pow( pow(r[0]*a_State.m_dx - 1.5, 2.0) + pow(r[1]*a_State.m_dx -1.5, 2.0), 2.0 )){

          sumVort +=  abs( *psi.data(r) - 1)*pow(a_State.m_dx, 2.0);

        } else{

           sumVort += abs( *psi.data(r))*pow(a_State.m_dx, 2.0);

        }


     }

     //Check conservation
     if( abs(a_State.wp - wg) > pow(10.0, -16.0)){

       cout << "conservation of vorticity not guaranteed after interpolation to grid in RK4 step at time " << a_time << " with error " << abs(wg-a_State.wp) << endl;

     }


   if( t_pt%step == 0){
     //Use Hockney to calculate streamfunction
     #ifdef PR_HDF5
      HDF5Handler h4;
      h4.writePatch(a_State.m_dx,psi, "vortg%i.hdf5", t_pt/step);
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


     //Error in velocity on grid
     for( auto it = B2.begin(); it != B2.end(); it++){
       r = it.operator*();

       z = pow( pow( (r[0]*a_State.m_dx - 1.5), 2.0) + pow( (r[1]*a_State.m_dx - 1.5), 2.0), 0.5 );

     if( z <= 1){

       sumUg += (abs(*Ug.data(r) - 0.5*(r[1]*a_State.m_dx - 1.5))+abs(*Vg.data(r) + 0.5*(r[0]*a_State.m_dx-1.5)) )*pow(a_State.m_dx, 2.0);

      } else{

        sumUg += (abs(*Ug.data(r) - 0.5/pow(z, 2.0)*(r[1]*a_State.m_dx - 1.5))+abs(*Vg.data(r) + 0.5/pow(z, 2.0)*(r[0]*a_State.m_dx-1.5)) )*pow(a_State.m_dx, 2.0);

      }

    }


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
     interpol.W22p( Ug, u, Y, a_State.hp, a_State.m_dx, Np);
     interpol.W22p( Vg, v, Y, a_State.hp, a_State.m_dx, Np);
     

     //velocity of particles error
     for(int i = 0; i < Np; i++){

         sumU +=  (abs(u[i] - 0.5*(a_State.X.at(i).x[1] - 1.5)) +abs(v[i] + 0.5*(a_State.X.at(i).x[0] - 1.5)))*pow(a_State.hp,2.0);
     }

     a_State.sumU = sumU;
     a_State.sumUg = sumUg;
     a_State.sumVort = sumVort;


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
    double strength, wc = 0;
    double k = 0;

    //Deposit on grid
    interpol.W22(strength_remap,X, hp, hp, Np);

    X.clear();
    particle p;

    // Generate new particle vector
    for(auto it = Bp.begin(); it != Bp.end(); it++){
      
      r = *it;

      strength = *strength_remap.data(r);
   
   
      if( strength > pow(10.0, -16.0 )){
        p.x[0] = r[0]*hp; p.x[1] = r[1]*hp;
	p.strength = strength;
        k++;     

	X.push_back(p);
      }

    }

    cout <<"Pre-remap # particles " <<  Np << " after remap # " << k << endl; 

    //Check conservation
    for( auto it = Bp.begin(); it != Bp.end(); it++){
       r = it.operator*();

       wc += *strength_remap.data(r)*pow(hp, 2.0);
    }

     if( abs(wp - wc) > pow(10.0, -17.0)){

       cout << "conservation of vorticity not guaranteed in remap " << abs(wc-wp) << endl;

     }

   
}
