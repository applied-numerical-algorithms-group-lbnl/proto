#ifndef _ADVECTION_H_
#define _ADVECTION_H_
#include "Proto.H"
#define NUMCOMPS 1

class Advection
{
public:
  inline void operator()(BoxData<double,NUMCOMPS>& a_UStage,
                         BoxData<double,NUMCOMPS>& a_kStage,
                         array<BoxData<double,NUMCOMPS>,DIM>& a_fluxes,
                         const double& a_dx)
  {
    // Advection operator for testing.
    array<double,DIM> uvec;
    uvec.fill(1.0);
    a_kStage.setVal(0.);
    double dxI = 1.0/a_dx;
    for (int dir = 0; dir < DIM; dir++)
      {
        Stencil<double> interp;
        if (uvec[dir] > 0)
          {
            interp = Stencil<double>::CellToFace(dir,Side::Lo,5);
          }
        else
          {
            interp = Stencil<double>::CellToFace(dir,Side::Hi,5)*((1.0)*Shift::Basis(dir,-1));
          }
        auto divergence = Stencil<double>::FluxDivergence(dir);
        a_fluxes[dir] |= interp(a_UStage,uvec[dir]);
        a_kStage -= divergence(a_fluxes[dir],dxI);
      }
  };
  static inline Point ghostSize(){return Point::Ones(3);};
};
  
#endif             
