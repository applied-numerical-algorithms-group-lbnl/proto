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
               const double a_dx, const double domain, const double w);
  double m_dx;
  double m_dt;
  double hp;
  double L;
  double wp;
  double sumU;
  double sumUg;
  double sumVort;
  vector<particle> X;
  void increment(const DX& a_DX);
  void getVelocity(BoxData<double>& vort,vector<double>& errorVg, vector<double>& errorpsi, vector<double>& u, vector<double>& v, const int Np,int step, State& a_State );
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
