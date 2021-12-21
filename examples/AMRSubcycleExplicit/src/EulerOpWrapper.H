#include "EulerOp.H"

using namespace EulerOp;
class EulerOpWrapper
{
public:
  inline void operator()(BoxData<double,NUMCOMPS>& a_U,
                         BoxData<double,NUMCOMPS>& a_divF,
                         array<BoxData<double,NUMCOMPS>,DIM>& a_fluxes,
                         BoxData<double,NUMCOMPS>& a_U0,
                         const double& a_dt,
                         const double& a_dx)
  {
    // Dummy operator for testing.
    a_divF.setVal(0.);
    for (int dir = 0; dir < DIM; dir++)
      {
        a_fluxes[dir].setVal(1.0);
      }
  };
  inline Point ghostSize(){return Point::Ones(2);};
};
  
             
