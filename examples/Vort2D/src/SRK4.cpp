#include "SRK4.H"
#include <iostream>
#include <vector>
#include "particle.hpp"
#include <cmath>
#include "interp.H"

State::State(const vector<particle> a,  const double a_dt,const double a_hp,  const double a_dx, const double domain)
{

        m_dt = a_dt;
	m_dx = a_dx;
        hp = a_hp;
	particle p;
	L = domain;

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

void F::operator()( DX& a_DX, double a_time, double a_dt, State& a_State){

     particle p;
          vector<particle> Y, vel;
     vector<double> u,v;
     vector<particle> Pg;
     int n = a_State.L/a_State.m_dx; //number of grid points - 1
     Box Bg(Point::Zeros(), Point(n,n) );   
     BoxData<double> psi(Bg);
     interp interpol;
//     int DIM = 2;
     int  Np = 0;
     int M = log2(n+1);

     Hockney hockney(a_State.m_dx, M);

     Point r;
     Stencil<double> m_derivative[DIM];

     for(unsigned int i = 0; i < a_DX.X2.size(); i++)
     {
          p.x[0] = a_DX.X2.at(i).x[0] + a_State.X.at(i).x[0];
	  p.x[1] = a_DX.X2.at(i).x[1] + a_State.X.at(i).x[1];
          p.strength = a_State.X.at(i).strength;
           
          Np++;
	  Y.push_back(p);

//	  cout << p.strength << endl;

     }

     //Interpolate to grid 
  //   cout << "interp begin" << endl;
     interpol.W22(psi, Y, a_State.m_dx, a_State.hp, Np); 
  //    cout << "interp end 1" << endl;

    //  for( auto it = Bg.begin(); it != Bg.end(); it++){
    //          r = it.operator*();
    //	      cout << psi(r) << endl;

    //  }
     //Use Hockney to calculate streamfunction

     hockney.convolve(psi);

//  for( auto it = Bg.begin(); it != Bg.end(); it++){
//              r = it.operator*();
//              cout << psi(r) << endl;

//      }


     //Determine velocity on the grid

     for( int dir = 0; dir < DIM; dir++){
        m_derivative[dir] = Stencil<double>::Derivative(1,dir,2);
     }
     BoxData<double> Ug(Bg); Ug  = m_derivative[1](psi);
     BoxData<double> Vg(Bg); Vg = -1*m_derivative[0](psi);
      
  for( auto it = Bg.begin(); it != Bg.end(); it++){
              r = it.operator*();
              cout << Ug(r)<< " " << Vg(r) << endl;

      }
  




     // vector of grid positions

//     cout << "interp 2" << endl;
     interpol.W22p( Ug, u, Y, a_State.hp, a_State.m_dx, Np);
     interpol.W22p( Vg, v, Y, a_State.hp, a_State.m_dx, Np);
     cout << "interp 2 end" << endl;
     
     for( int i = 0;  i < Np; i++){

        p.x[0] = a_dt*u[i]; p.x[1] =a_dt*v[i];
	p.strength = 0;

	vel.push_back(p);
        cout << p.x[0] << " " << p.x[1] << endl;

     }

   

     a_DX.X2 = vel;

     

 }
