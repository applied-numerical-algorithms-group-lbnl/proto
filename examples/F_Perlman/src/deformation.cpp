#include "Proto.H"
#include <iostream>
#include <particle.hpp>
#include <vector>
#include <deform.hpp>
#include "deformation.H"
#include "SRK4.H"
#include "Proto_Timer.H"

using namespace Proto;
deformation::deformation(const vector<particle> y){
PR_TIMERS("main::RK4_f::deformation_init");
        deform d;

	for(int i = 0; i < y.size(); i++){

          d.F[0][0] = 1; d.F[1][1] = 1;
	  d.F[1][0] = 0; d.F[0][1] = 0;

	  F.push_back(d);

	  X.push_back(y.at(i));
	 

	}

}

void deformation::update_X( State& state){

PR_TIMERS("main::RK4_f::updateX");
	
	X = state.X;
        a_State.X = state.X;
	a_State.L = state.L;
	a_State.m_dx = state.m_dx;
	a_State.hp = state.hp;

}

void
DF::increment(double& a_weight,const DF& a_incr)
{

   PR_TIMERS("main::RK4_f::increment");
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

    PR_TIMERS("main::RK4_f::F*");
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
     Stencil<double> first_derivative[DIM];
     Stencil<double> mixed_derivative;
     long double z;
     long double sumpsi = 0;
     double sumUg = 0;
      
     psi.setToZero();
     Hockney hockney(a_State.m_dx, M);

     interpol.W44(psi, a_State.X, a_State.m_dx, a_State.hp, Np);

     hockney.convolve(psi);
 
     for( int dir = 0; dir < DIM; dir++){
       m_derivative[dir] = Stencil<double>::LaplacianFace(dir,2);
       first_derivative[dir] = Stencil<double>::Derivative(1, dir, 2);
     }

   
     mixed_derivative = first_derivative[0]*first_derivative[1];

     dudxg.setToZero(); dudyg.setToZero(); 
     dvdxg.setToZero(); dvdyg.setToZero();

     //du_jdx_i calculations
     dudyg = m_derivative[0](psi, 1.0/pow(a_State.m_dx, 2.0)); dudxg = mixed_derivative(psi, 1.0/pow(a_State.m_dx, 2.0));
     dvdxg = m_derivative[1](psi, 1.0/pow(a_State.m_dx, 2.0)); 
     dvdxg *= -1;
     dudxg.copyTo(dvdyg);
     dvdyg *= -1;
   
    #ifdef PR_HDF5
       HDF5Handler h5;
       h5.writePatch(a_State.m_dx,dudxg, "dudxg%i.hdf5",0);
       h5.writePatch(a_State.m_dx,dudyg, "dudyg%i.hdf5",0);
       h5.writePatch(a_State.m_dx,dvdxg, "dvdxg%i.hdf5",0);
       h5.writePatch(a_State.m_dx,dvdyg, "dvdyg%i.hdf5",0);

    #endif

    //Interpolate from grid to particles   
    interpol.W44p( dudxg, dudx, Y, a_State.hp, a_State.m_dx, Np);
    interpol.W44p( dvdxg, dvdx, Y, a_State.hp, a_State.m_dx, Np);

    interpol.W44p( dudyg, dudy, Y, a_State.hp, a_State.m_dx, Np);
    interpol.W44p( dvdyg, dvdy, Y, a_State.hp, a_State.m_dx, Np);
     
}
  
void RHS::operator()( DF& a_DF, double a_time, double a_dt, deformation& state){
     
     deform d;
     int n = state.a_State.L/state.a_State.hp; //number of grid points - 1
     vector<particle> Y;
     vector<double> dudx, dvdx, dudy, dvdy;
     particle p;
     int Np = state.F.size();

     
     { //Get dudx 
     PR_TIMERS("main::RK4_f::getv");
     getVelocity(dudx, dvdy, dudy, dvdx, Np, state.a_State, state.X);
     }


      //Update what k should be
      for(unsigned int i = 0; i < state.F.size(); i++){

	PR_TIMERS("main::RK4_f::updateF");
        d.F[0][0] = (dudx[i]*(state.F.at(i).F[0][0] + a_DF.F2.at(i).F[0][0]) + dudy[i]*(state.F.at(i).F[1][0] + a_DF.F2.at(i).F[1][0]))*a_dt;
        d.F[0][1] = (dudx[i]*(state.F.at(i).F[0][1] + a_DF.F2.at(i).F[0][1]) + dudy[i]*(state.F.at(i).F[1][1] + a_DF.F2.at(i).F[1][1]))*a_dt;
	d.F[1][0] = (dvdx[i]*(state.F.at(i).F[0][0] + a_DF.F2.at(i).F[0][0]) + dvdy[i]*(state.F.at(i).F[1][0] + a_DF.F2.at(i).F[1][0]))*a_dt;
	d.F[1][1] = (dvdx[i]*(state.F.at(i).F[0][1] + a_DF.F2.at(i).F[0][1]) + dvdy[i]*(state.F.at(i).F[1][1] + a_DF.F2.at(i).F[1][1]))*a_dt;

        a_DF.F2.at(i).F[0][0] = d.F[0][0]; 
        a_DF.F2.at(i).F[0][1] = d.F[0][1];
        a_DF.F2.at(i).F[1][0] = d.F[1][0];
        a_DF.F2.at(i).F[1][1] = d.F[1][1];

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
 
 void deformation::QR(vector<deform>& Q, vector<deform>& R, vector<double>& angleQ, vector<double>& eigenR ){

	 double e1[2], e2[2], u1[2], u2[2],v1[2], v2[2];
	 deform Q_m, R_m, QQT;
	 double uv;
         double v1_norm, u2_norm;

	 Q.clear(); R.clear(); angleQ.clear(); eigenR.clear(); 


         //Gram Schmidt
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

	   //Construct R
	   R_m.F[0][0] = v1_norm; R_m.F[0][1] = uv;
	   R_m.F[1][0] = 0; R_m.F[1][1] = u2_norm;

	   //Construct Q
	   Q_m.F[0][0] = e1[0]; Q_m.F[1][0] = e1[1];
	   Q_m.F[0][1] = e2[0]; Q_m.F[1][1] = e2[1];
 
	  // Q & R determinant 
	  //f( i == 0){
          //Q Determinant
          //if( abs( (Q_m.F[0][0]*Q_m.F[1][1] - Q_m.F[0][1]*Q_m.F[1][0]) -1 ) > pow(10.0, -6.0) ){
          //cout << "det(Q) =/= 1 " << abs(Q_m.F[0][0]*Q_m.F[1][1] - Q_m.F[0][1]*Q_m.F[1][0]-1)  << endl;
          //
          //}
          //if( abs( (R_m.F[0][0]*R_m.F[1][1] - R_m.F[0][1]*R_m.F[1][0]) -1 ) > pow(10.0, -6.0) ){
          //cout << "det(R) =/= 1 " << abs(R_m.F[0][0]*R_m.F[1][1] - R_m.F[0][1]*R_m.F[1][0] - 1)  << endl;
          //cout << R_m.F[0][0] << " " << R_m.F[0][1] << endl;
          //cout << R_m.F[1][0] << " " << R_m.F[1][1] << endl;
          //}
          //}

	   // QQT
	   QQT.F[0][0] = pow( Q_m.F[0][0], 2.0) + pow( Q_m.F[0][1], 2.0);
	   QQT.F[0][1] = Q_m.F[0][0]*Q_m.F[1][0] + Q_m.F[1][1]*Q_m.F[0][1];
	   QQT.F[1][0] = Q_m.F[0][0]*Q_m.F[1][0] + Q_m.F[1][1]*Q_m.F[0][1];
	   QQT.F[1][1] = pow( Q_m.F[1][1], 2.0) + pow( Q_m.F[1][0], 2.0);


           //Check QQT = 1	   
	   if( (abs(QQT.F[0][0] - 1) > pow( 10.0, -7.0) ) || (abs(QQT.F[1][1] - 1) > pow( 10.0, -7.0) ) || (abs(QQT.F[0][1] ) > pow( 10.0, -7.0) ) ){

	      cout << "QQT =/= I " << endl;
              cout << QQT.F[0][0] << " " << QQT.F[0][1] << endl;
	      cout << QQT.F[1][0] << " " << QQT.F[1][1] << endl;
	   }
				  
	   Q.push_back(Q_m);
	   R.push_back(R_m);


           if( v1_norm >= 1){

	     eigenR.push_back(v1_norm);

           }else{
             
	     eigenR.push_back(u2_norm);

	   }

	   angleQ.push_back( acos(Q_m.F[0][0]) );

	 }


   }

void deformation::remap(State& a_state){


   F.clear();
   X.clear();
   deform d;

   for(int i = 0; i < a_state.X.size(); i++){

      d.F[0][0] = 1; d.F[1][1] = 1;
      d.F[1][0] = 0; d.F[0][1] = 0;

      F.push_back(d);
      X.push_back(a_state.X.at(i));


   }

}
