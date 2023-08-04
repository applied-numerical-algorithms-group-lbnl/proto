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
        diffVdt = 0;
        diffVort = 0;	

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
     double c = 0.75 ;
     long double sumpsi = 0;
     double sumUg = 0, sumUg1 = 0;
     vector<double> vortpgp;
     //Interpolate from grid to particles
     interpol.W44(vort, a_State.X, a_State.m_dx, a_State.hp, Np);
     double h2 = pow(a_State.m_dx, 2.0); 
     int px, py;

     Point i_closest;
     psi.setToZero();
     a_State.radi.clear();
     a_State.corrVort_error.clear();
     a_State.corrRelVort.clear();
     a_State.vortpgp.clear();
     a_State.vortpgpdt.clear();
     a_State.vortpgp_error.clear();
     a_State.vortpgpdt_error.clear();


     //Calculate relative and exact error in vorticity
     for(int i = 0; i < Np; i++){

       px = static_cast<int> ( round((a_State.X.at(i).x[0]+c) / a_State.m_dx) );
       py = static_cast<int> ( round((a_State.X.at(i).x[1]+c) / a_State.m_dx) );

       i_closest = Point({px, py} );

     
       a_State.corrVort_error.push_back( abs(*vort.data(i_closest)-a_State.X.at(i).strength) );
       a_State.corrRelVort.push_back( abs(*vort.data(i_closest)-a_State.X.at(i).strength) / abs(a_State.X.at(i).strength) ); 

     }


     vort.copyTo(psi, B, B); hockney.convolve(psi); 
    
     //Calculate velocity from Hockney results
      for( int dir = 0; dir < DIM; dir++){
        m_derivative[dir] = Stencil<double>::Derivative(1,dir,4);
     }
     
     BoxData<double> Ug(B); Ug.setToZero();
     Ug  = m_derivative[1](psi, 1.0/a_State.m_dx);
     BoxData<double> Vg(B); Vg.setToZero();
     Vg = m_derivative[0](psi, 1.0/a_State.m_dx);
     Vg *= -1;
     Point r;


     //Interpolate to get particle velocity
     interpol.W44p( Ug, u, a_State.X, a_State.hp, a_State.m_dx, Np);
     interpol.W44p( Vg, v, a_State.X, a_State.hp, a_State.m_dx, Np);
     interpol.W44p( vort, vortpgp, a_State.X, a_State.hp, a_State.m_dx, Np); 

     a_State.sumVortpgp = 0;
     a_State.sumVortpgpdt = 0;
     for(int i = 0; i < Np; i++){


         a_State.vortpgp_error.push_back( abs(vortpgp.at(i) - a_State.X.at(i).strength ) );
	 a_State.vortpgpdt_error.push_back( 0 );
         a_State.vortpgpdt.push_back(0);
	 a_State.vortpgp.push_back(vortpgp.at(i));

	 a_State.sumVortpgp += pow( vortpgp.at(i) - a_State.X.at(i).strength, 2.0)*pow(a_State.hp, 2.0);

     }

         a_State.sumVortpgp = pow(a_State.sumVortpgp, 0.5);


}


void F::operator()( DX& a_DX, double a_time, double a_dt, State& a_State){
    
     particle p;
     particle lam;
     vector<particle> Y, vel;
     vector<double> u,v, vortpgp;
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
     BoxData<double> vort(Bg);
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
     double c = 0.75 ;
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
 

     psi.copyTo(vort, Bg, Bg);

     //Calculate error in vorticity field
     for(int i = 0; i < Np; i++){
       px = static_cast<int> ( round((a_State.X.at(i).x[0]+c) / a_State.m_dx) );
       py = static_cast<int> ( round((a_State.X.at(i).x[1]+c) / a_State.m_dx) );

       i_closest = Point({px, py} );

       a_State.corrVort_error.push_back( abs(*psi.data(i_closest)-a_State.X.at(i).strength) );

       a_State.corrRelVort.push_back(abs(*psi.data(i_closest)-a_State.X.at(i).strength)/abs(a_State.X.at(i).strength) );

     }
     array<double,DIM> temp;
     BoxData<double> lambda_error(Bg); lambda_error.setToZero();
     //Calculate voriticty summation over grid and error
     for( auto it = Bg.begin(); it != Bg.end(); it++){
       r = *it;


       wg += *psi.data(r)*pow(a_State.m_dx, 2.0);

       lamb_g(r) = (*psi.data(r))*(*lambda_g.data(r));
       lambda_error(r) += abs(abs(*lambda_g.data(r) -1)*(*psi.data(r)));

     }

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
     
    string Ugfile = "Ug_" + to_string(t_pt/step) + ".txt";
    ofstream f1(Ugfile);
    string Vgfile = "Vg_" + to_string(t_pt/step) + ".txt";
    ofstream f2(Vgfile);


    if( f1.is_open() ) {

      //format in scientific notation
       f1.setf(std::ios::scientific, std::ios::floatfield);
       f1.precision(8);


      for( auto it = B2.begin(); it != B2.end(); it++){
       r = *it;


       f1 << *Ug.data(r) << endl;


     }


   }

   if( f2.is_open() ) {

      //format in scientific notation
       f2.setf(std::ios::scientific, std::ios::floatfield);
       f2.precision(8);


      for( auto it = B2.begin(); it != B2.end(); it++){
       r = *it;


       f2 << *Vg.data(r) << endl;


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
 


     { //Interpolate from grid to particles

     PR_TIMERS("main::RK4::interp_grid");
     interpol.W44p( Ug, u, Y, a_State.hp, a_State.m_dx, Np);

     interpol.W44p( Vg, v, Y, a_State.hp, a_State.m_dx, Np);
     interpol.W44p( vort, vortpgp, a_State.X, a_State.hp, a_State.m_dx, Np);

     }

     a_State.sumVortpgp = 0;
     a_State.vortpgp_error.clear();
     a_State.vortpgp = vortpgp;

     if( Np > a_State.vortpgpdt_error.size() ){
         
        a_State.vortpgpdt.clear();
	a_State.vortpgpdt_error.clear();

     }
     for(int i = 0; i < Np; i++){


         a_State.vortpgp_error.push_back( abs(vortpgp.at(i) - a_State.X.at(i).strength ) );
         a_State.sumVortpgp += pow( vortpgp.at(i) - a_State.X.at(i).strength, 2.0)*pow(a_State.hp, 2.0);

	 if( Np > a_State.vortpgpdt.size()){
  
             a_State.vortpgpdt_error.push_back( a_State.vortpgp_error.at(i)*a_dt);
             a_State.vortpgpdt.push_back(a_State.vortpgp.at(i)*a_dt);


         }else{

             a_State.vortpgpdt_error.at(i) = a_State.vortpgpdt_error.at(i) + a_State.vortpgp_error.at(i)*a_dt;
             a_State.vortpgpdt.at(i) = a_State.vortpgpdt.at(i) + a_State.vortpgp.at(i)*a_dt;

	   

	 }
     }


     a_State.sumVortpgp = pow(a_State.sumVortpgp, 0.5);
     a_State.sumVortpgpdt = a_State.sumVortpgpdt + a_dt*a_State.sumVortpgp;


     //generate updated info for k
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
    double c = 0.75 ;
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
  
      if( abs(strength) > pow(10, -7.0)){
	      
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

void State::rate(State& second_State, int step ){

     interp interpol;
     int n = L/m_dx; //number of grid points - 1
     int px, py,px2,py2;
     int M = log2(n+1);
     Box Bg(Point::Zeros(), Point(n,n) );
     Box B2( Point({2,2}), Point({n-2, n-2}) );
     BoxData<double> psi(Bg);
     BoxData<double> psi_2nd(Bg);
     BoxData<double> psiDiff(Bg);
     BoxData<double> Ug(Bg), Vg(Bg),Ug_2nd(Bg), Vg_2nd(Bg);
     BoxData<double> Ug_diff(B2), Vg_diff(B2), Vel_diff(B2);
     Point r;
     int c = 0.75 ;
     Stencil<double> m_derivative[DIM];
     int Np = X.size();
     int Np2 = second_State.X.size();
     psi.setToZero();
     psi_2nd.setToZero();
     psiDiff.setToZero();
     interpol.W44(psi, X, m_dx, hp, Np);
     interpol.W44(psi_2nd,second_State.X, m_dx, second_State.hp, Np2);
     
     diffV1   = 0;
     diffV2   = 0;
     diffVort = 0;
     diffUg   = 0;
     diffVg   = 0;
     diffV    = 0;
     Hockney hockney(m_dx, M);
     for( auto it =  Bg.begin(); it != Bg.end(); it++){

	 r = *it;
	 
         psiDiff(r) += abs(psi(r) - psi_2nd(r) );
         diffVort += pow((psi(r) - psi_2nd(r))*m_dx, 2.0);


     }

     diffVort = pow(diffVort, 0.5);
     diffVortdt += m_dt*diffVort;
     #ifdef PR_HDF5
     HDF5Handler h5;
     h5.writePatch(m_dx,psi, "vort_rate_%i.hdf5", step);
     h5.writePatch(m_dx,psi_2nd, "vort_rate_2_diff_%i.hdf5", step);
     h5.writePatch(m_dx,psiDiff, "Vort_diff_%i.hdf5", step);
     #endif
 
     //Determine velocity on the grid
     for( int dir = 0; dir < DIM; dir++){
        m_derivative[dir] = Stencil<double>::Derivative(1,dir,4);
     }

     hockney.convolve(psi);
     hockney.convolve(psi_2nd);
     Ug.setToZero(); Ug_2nd.setToZero();
     Ug  = m_derivative[1](psi, 1.0/m_dx);
     Ug_2nd  = m_derivative[1](psi_2nd, 1.0/m_dx);
     Vg.setToZero(); Vg_2nd.setToZero();
     Vg = m_derivative[0](psi, 1.0/m_dx );
     Vg_2nd = m_derivative[0](psi_2nd, 1.0/m_dx );
     Vg *= -1; Vg_2nd *= -1;

     Vg_diff.setToZero(); Ug_diff.setToZero();
     Vel_diff.setToZero();
     for( auto it = B2.begin(); it != B2.end(); it++){

         r = *it;
         Vg_diff(r) += abs(Vg(r) - Vg_2nd(r) );
	 Ug_diff(r) += abs(Ug(r) - Ug_2nd(r) );
	 Vel_diff(r) += pow( pow((Ug(r) - Ug_2nd(r)), 2.0) + pow((Vg(r) - Vg_2nd(r)), 2.0), 0.5);
         diffV1 += (  abs(Vg(r) - Vg_2nd(r) ) + abs(Ug(r) - Ug_2nd(r) ))*pow(m_dx, 2.0);
	 diffV2 += pow((Ug(r) - Ug_2nd(r)), 2.0) + pow((Vg(r) - Vg_2nd(r))*m_dx, 2.0); 
         

     }

     string Vfile = "VelDiff_" + to_string(m_dt/step) + ".txt";
     ofstream f1(Vfile);


     if( f1.is_open() ) {

      //format in scientific notation
       f1.setf(std::ios::scientific, std::ios::floatfield);
       f1.precision(8);


        for( auto it = B2.begin(); it != B2.end(); it++){
         r = *it;


         f1 << *Vel_diff.data(r) << endl;


       }


     }

     diffVg = Vg_diff.max();
     diffUg = Ug_diff.max();
     diffV  = Vel_diff.max();
     diffV2 = pow(diffV2, 0.5);
     diffVdt += m_dt*Vel_diff.max();
     #ifdef PR_HDF5
     h5.writePatch(m_dx,Ug_diff, "Ug_diff_%i.hdf5", step);
     h5.writePatch(m_dx,Vg_diff, "Vg_diff_%i.hdf5", step);
     h5.writePatch(m_dx,psiDiff, "Vort_diff_%i.hdf5", step);
     h5.writePatch(m_dx,Vel_diff, "Vel_diff_%i.hdf5", step);
     h5.writePatch(m_dx,Ug, "Ug_1st_%i.hdf5", step);
     h5.writePatch(m_dx,Vg, "Vg_1st_%i.hdf5", step);
     h5.writePatch(m_dx,Ug_2nd, "Ug_2nd_%i.hdf5", step);
     h5.writePatch(m_dx,Vg_2nd, "Vg_2nd_%i.hdf5", step);
     #endif





}


