#ifndef _interp_H_
#define _interp_H_

#include "Proto.H"
#include <iostream>
#include <particle.hpp>
#include <vector>

using namespace Proto;



class interp
{
public:
 void W22(BoxData<double>& A, const vector<particle> p,
               const double h_to,
               const double h_from,
	       const int Np);
 void W22p(BoxData<double>& A,  vector<double>& p,
	       const vector<particle> X, 
	       const double h_to,
	       const double h_from,
	       const int Np );

 void W44(BoxData<double>& A, const vector<particle> p,
               const double h_to,
               const double h_from,
               const int Np);
 void W44p(BoxData<double>& A,  vector<double>& p,
               const vector<particle> X,
               const double h_to,
               const double h_from,
               const int Np );

};

#endif
