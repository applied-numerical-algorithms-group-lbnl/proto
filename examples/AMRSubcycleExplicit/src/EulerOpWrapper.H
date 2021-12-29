#ifndef _ADVECTION_H_
#define _ADVECTION_H_
#include "Proto.H"
#define NUMCOMPS DIM+2

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
        interp = Stencil<double>::CellToFace(dir,Side::Lo,5);
#if 1
        interp = (7./12.)*Shift(Point::Basis(dir,-1))
          + (7./12.)*Shift(Point::Zeros())
          + (-1./12.)*Shift(Point::Basis(dir,-2))
          + (-1./12.)*Shift(Point::Basis(dir,1));
#endif
#if 0
        interp =
          (2./60.)*Shift::Basis(dir,-3) +
          (-13./60.)*Shift::Basis(dir,-2) +
          (47./60.)*Shift::Basis(dir,-1) +
          (27./60.)*Shift::Zeros() +
          (-3./60.)*Shift::Basis(dir,1);
#endif
        auto divergence = Stencil<double>::FluxDivergence(dir);
        a_fluxes[dir] |= interp(a_UStage,uvec[dir]);
        a_kStage -= divergence(a_fluxes[dir],dxI);
      }
  };
  static inline Point ghostSize(){return Point::Ones(3);};
};
  
#endif             
