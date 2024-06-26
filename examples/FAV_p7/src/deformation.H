#ifndef _deformation_H_
#define _deformation_H_

#include "Proto.H"
#include <iostream>
#include <particle.hpp>
#include <vector>
#include <deform.hpp>
#include "SRK4.H"
using namespace Proto;


class DF;
class deformation
{
public:
   deformation(){};
  ~deformation(){};
   deformation(const vector<particle> x);

   vector<deform> F;
   vector<particle> X;
   State a_State;
   void update_X( State& a_State);
   
   void increment(const DF &a_DF);
   void remap(State& a_state); 
   void QR(vector<deform>& Q, vector<deform>& R, vector<double>& angleQ, vector<double>& eigenR, vector<double>& detF );
};

class DF
{
public:
  DF(){};
  ~DF(){};
  vector<deform> F2;
  void increment(double& a_weight,const DF& a_incr);
  void init(deformation& state);
  void operator*=(const double& a_weight);
};


class RHS
{
public:
      RHS(){};
      ~RHS(){};
      void operator()( DF& a_DF, double a_time, double a_dt, deformation& F);
      void getVelocity(vector<double>& dudx, vector<double>& dvdy, vector<double>& dudy, vector<double>& dvdx, const int Np, State& a_State, const vector<particle> Y );
};

#endif
