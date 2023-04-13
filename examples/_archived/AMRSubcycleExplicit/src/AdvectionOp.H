#ifndef _ADVECTION_H_
#define _ADVECTION_H_
#include "Proto.H"
#define NUMCOMPS 1
PROTO_KERNEL_START
void
thresholdF(
           Var<short,1>& a_tags,
           Var<double,NUMCOMPS>& a_U)
{
  double thresh = 1.0001;
  if (a_U(0) > thresh)
    {a_tags(0) = 1;}
  else
    {a_tags(0) = 0;};
};
PROTO_KERNEL_END(thresholdF,threshold);
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
        BoxData<double,NUMCOMPS> temp = divergence(a_fluxes[dir],dxI);
        a_kStage -= temp;
        PROTO_ASSERT(temp.box().contains(a_kStage.box()),"Too few ghost cells for divergence");
        //a_kStage -= divergence(a_fluxes[dir],dxI);
      }
  };
  static inline Point ghostSize(){return Point::Ones(3);}
  static inline int order(){return 5;}
  static inline void generateTags(LevelTagData& a_tags,
                                  LevelBoxData<double,NUMCOMPS>& a_U,
                                  double a_dx,
                                  double a_time)
  {
    for (auto dit = a_tags.begin(); *dit != dit.end();++dit)
      {
        forallInPlace(threshold,a_tags[*dit],a_U[*dit]);
      }
  }
};  
#endif             
