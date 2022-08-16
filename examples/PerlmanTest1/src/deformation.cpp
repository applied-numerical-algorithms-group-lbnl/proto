#include "Proto.H"
#include <iostream>
#include <particle.hpp>
#include <vector>
#include <deform.hpp>
#include "deformation.H"
#include "SRK4.H"
using namespace Proto;
deformation::deformation(const vector<particle> y){

        deform d;
	for(int i = 0; i < y.size(); i++){

          d.F[0][0] = 1; d.F[1][1] = 1;
	  d.F[1][0] = 0; d.F[0][1] = 0;

	  F.push_back(d);

	  X.push_back(y.at(i));
	 

	}

}

void deformation::update_X( State& state){

	X = state.X;
	a_State = state;

}

void
DF::increment(double& a_weight,const DF& a_incr)
{
   for(unsigned int i = 0; i < F2.size(); i++){
     F2.at(i).F[0][0] += a_weight*a_incr.F2.at(i).F[0][0];
     F2.at(i).F[0][1] += a_weight*a_incr.F2.at(i).F[0][1];
     F2.at(i).F[1][0] += a_weight*a_incr.F2.at(i).F[1][0];
     F2.at(i).F[1][1] += a_weight*a_incr.F2.at(i).F[1][1];

   }

   

};

void
DF::init(deformation& state)
{
    deform d;

    F2 = state.F;

    for(unsigned int i = 0; i < state.F.size(); i++){

            F2.at(i).F[0][0] = 0;
            F2.at(i).F[0][1] = 0;
            F2.at(i).F[1][0] = 0;
            F2.at(i).F[1][1] = 0;


    }

   

};

void
DF::operator*=(const double& a_weight)
{

    for(unsigned int i = 0; i < F2.size(); i++){
        F2.at(i).F[0][0] *= a_weight;
        F2.at(i).F[0][1] *= a_weight;
	F2.at(i).F[1][0] *= a_weight;
        F2.at(i).F[1][1] *= a_weight;

    }

   
}


void RHS::getVelocity(vector<double>& dudx, vector<double>& dvdy, vector<double>& dudy, vector<double>& dvdx, const int Np, State& a_State, const vector<particle> Y ){

     int n = a_State.L/a_State.m_dx; //number of grid points - 1
     Box Bg(Point::Zeros(), Point(n,n) );
     Box B2(Point(1,1), Point(n-1, n-1) );
     int np = a_State.L/a_State.hp;
     BoxData<double> psi(Bg);
     BoxData<double> dudxg(Bg), dudyg(Bg), dvdxg(Bg), dvdyg(Bg);
     interp interpol;
     int M = log2(n+1);
     Stencil<double> m_derivative[DIM];
     Hockney hockney(a_State.m_dx, M);

     interpol.W22(psi, a_State.X, a_State.m_dx, a_State.hp, Np);

     hockney.convolve(psi);
 
     for( int dir = 0; dir < DIM; dir++){
        m_derivative[dir] = Stencil<double>::Derivative(1,dir,2);
     }

     BoxData<double> Ug(Bg);
     Ug.setToZero();
     Ug  = m_derivative[1](psi);
     BoxData<double> Vg(Bg);
     Vg.setToZero();
     Vg = m_derivative[0](psi);
     Vg *= -1;

     dudxg.setToZero(); dudyg.setToZero(); 
     dvdxg.setToZero(); dvdyg.setToZero();

     dudxg = m_derivative[0](Ug); dudyg = m_derivative[1](Ug);
     dvdyg = m_derivative[1](Vg); dvdxg = m_derivative[0](Vg);

     interpol.W22p( dudxg, dudx, Y, a_State.hp, a_State.m_dx, Np);
     interpol.W22p( dvdxg, dvdx, Y, a_State.hp, a_State.m_dx, Np);

     interpol.W22p( dudyg, dudy, Y, a_State.hp, a_State.m_dx, Np);
     interpol.W22p( dvdyg, dvdy, Y, a_State.hp, a_State.m_dx, Np);

     
}
  
void RHS::operator()( DF& a_DF, double a_time, double a_dt, deformation& state){

     int n = state.a_State.L/state.a_State.hp; //number of grid points - 1
     vector<particle> Y;
     vector<double> dudx, dvdx, dudy, dvdy;
     particle p;
     int Np = state.F.size();
         
      getVelocity(dudx, dvdy, dudy, dvdx, Np, state.a_State, state.X);
      

      for(int i = 0; i < state.F.size(); i++){

        a_DF.F2.at(i).F[0][0] += (dudx[i]*(state.F.at(i).F[0][0]+a_DF.F2.at(i).F[0][0]) + dudy[i]*(state.F.at(i).F[1][0] + a_DF.F2.at(i).F[1][0] ))*a_dt;
        a_DF.F2.at(i).F[0][1] += (dudx[i]*(state.F.at(i).F[0][1]+a_DF.F2.at(i).F[0][1]) + dudy[i]*(state.F.at(i).F[1][1] + a_DF.F2.at(i).F[1][1]))*a_dt;
	a_DF.F2.at(i).F[1][0] += (dvdx[i]*(state.F.at(i).F[0][0] + a_DF.F2.at(i).F[0][0]) + dvdy[i]*(state.F.at(i).F[1][0] + a_DF.F2.at(i).F[1][0]))*a_dt;
	a_DF.F2.at(i).F[1][1] += (dvdx[i]*(state.F.at(i).F[0][1] + a_DF.F2.at(i).F[0][1]) + dvdy[i]*(state.F.at(i).F[1][1] + a_DF.F2.at(i).F[1][1]))*a_dt;

      }



}


void deformation::increment(const DF& a_DF){


	for(int i = 0; i < a_DF.F2.size(); i++){

		F.at(i).F[0][0] += a_DF.F2.at(i).F[0][0];
		F.at(i).F[0][1] += a_DF.F2.at(i).F[0][1];
		F.at(i).F[1][0] += a_DF.F2.at(i).F[1][0];
		F.at(i).F[1][1] += a_DF.F2.at(i).F[1][1];
	}



}
 
 void deformation::QR(vector<deform>& Q, vector<deform>& R ){

	 double e1[2], e2[2], u1[2], u2[2],v1[2], v2[2];
	 deform Q_m, R_m;
	 double uv;
         double v1_norm, u2_norm;
         Q.clear(); R.clear();
	 for( int i = 0; i < F.size(); i++){

           v1[0] = F.at(i).F[0][0];
	   v1[1] = F.at(i).F[1][0];
	   v2[0] = F.at(i).F[0][1];
           v2[1] = F.at(i).F[1][1];

	   v1_norm = pow( pow(v1[0], 2.0) + pow(v1[1], 2.0), 0.5);

	   e1[0] = v1[0]/v1_norm; e1[1] = v1[1]/v1_norm;

	   uv = v2[0]*e1[0] + v2[1]*e1[1];

	   u2[0] = v2[0] - uv*e1[0]; u2[1] = v2[1] - uv*e1[1];


	   u2_norm = pow( pow(u2[0], 2.0) + pow(u2[1], 2.0), 0.5);

	   e2[0] = u2[0]/u2_norm; e2[1] = u2[1]/u2_norm;

	   R_m.F[0][0] = v1_norm; R_m.F[0][1] = uv;
	   R_m.F[1][0] = 0; R_m.F[1][1] = u2_norm;

	   Q_m.F[0][0] = e1[0]; Q_m.F[1][0] = e1[1];
	   Q_m.F[0][1] = e2[0]; Q_m.F[1][1] = e2[1];

	   Q.push_back(Q_m);
	   R.push_back(R_m);



	 }


   }



