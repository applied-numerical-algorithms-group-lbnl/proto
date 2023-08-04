#ifndef _SRK4_H_
#define _SRK4_H_

#include "Proto.H"
#include <iostream>
#include <particle.hpp>
#include <vector>
#include "Hockney.H"
#include "interp.H"

using namespace Proto;



class DX;
class State
{
public:
  State(){};
  ~State(){};
   State(const vector<particle> a,
               const double a_dt,
	       const double a_hp,
               const double a_dx, const double domain, const double w, const vector<double> omega);
  double m_dx;
  double m_dt;
  double hp;
  double L;
  double wp;
  vector<double> radi;
  vector<double> lambda;
  double sumU, sumUm, sumU1;
  double sumUg, sumUg1, sumUgm;
  vector<double> corrVort_error;
  vector<double> corrRelVort;
  vector<double> vortpgp, vortpgp_error;
  vector<double> vortpgpdt, vortpgpdt_error;
  double sumVort, sumVort1, sumVortm;
  double sumVortpgpdt, sumVortpgp;
  double diffVort,diffUg,diffVg, diffV, diffVdt,diffVortdt;
  double diffV1, diffV2;
  vector<particle> X;
  BoxData<double> Vort_error;
  void increment(const DX& a_DX);
  void remap();
  void getVelocity(BoxData<double>& vort,vector<double>& errorVg, vector<double>& errorpsi, vector<double>& u, vector<double>& v, const int Np,int step, State& a_State );
  void rate(State& second_state, int step);

};

class DX
{
public:
  DX(){};
  ~DX(){};
  vector<particle> X2;
  void increment(double& a_weight,const DX& a_incr);
  void init(State& a_State);
  void operator*=(const double& a_weight);
};

class F
{
public:
      F(){};
      ~F(){};
      void operator()(
                 DX& a_DX,
                 double a_time,
                 double a_dt,
                 State& a_State);
};

#endif
