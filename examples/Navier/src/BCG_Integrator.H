#ifndef _BCG_Integrator_H_
#define _BCG_Integrator_H_

#include "Proto.H"


using namespace Proto;
using std::string;
///
/**
 */
class BCG_Integrator
{
public:
  static double s_dx;
  /// Default constructor.

  BCG_Integrator(const Box    & a_domain,
                 const double & a_dx,
                 const double & a_viscosity,
                 const int    & a_nghost)
  {
    m_domain    = a_domain;
    m_dx        = a_dx;
    m_viscosity = a_viscosity;
    m_nghost = a_nghost;
    defineStencils();
  }

  //divides cell centered velocity into vortical compoenent and a gradient of a scalar
  void
  ccProject(BoxData<double, DIM>& a_velocity,
            BoxData<double, DIM>& a_gradpres);

  ///
  /**
     Given velocity at n and gradient of pressure at n-1/2, advance both using BCG
   */
  void advanceSolution(BoxData<double, DIM>& a_velocity, 
                       BoxData<double, DIM>& a_gradpres, 
                       const double        & a_dt);

#if DIM==2
  ///
  void 
  computeVorticity(BoxData<double,   1>& vorticity,
                   BoxData<double, DIM>& velocity);
#else
  ///
  void 
  computeVorticity(BoxData<double, DIM>& vorticity,
                   BoxData<double, DIM>& velocity);
#endif

  ///public for testing
  void
  getUDotDelU(BoxData<double, DIM> & a_udelu,
              BoxData<double, DIM> & a_velocity,
              const double & a_dt);

  ///  Destructor.
  ~BCG_Integrator(){};


private:
  //weak construction forbidden
  BCG_Integrator();  


  ///fills ghost cells of velocity with periodic image
  void
  enforceBoundaryConditions(BoxData<double, DIM>& a_velocity);


  //divides face centered velocity into vortical compoenent and a gradient of a scalar
  void
  MACProject(BoxData<double,   1> a_velocity[DIM],
             BoxData<double, DIM> a_gradpres[DIM]);


  //computes nu*lapl(velcomp)
  void
  getSourceTerm(BoxData<double, 1> & a_source,
                BoxData<double, 1> & a_velcomp);
  


  void
  averageVelocityToFaces(BoxData<double,   1>   a_velface[DIM],
                         BoxData<double, DIM> & a_velcell);


  void
  MACDivergence(BoxData<double, 1> & a_divergence,
                BoxData<double, 1>   a_velface[DIM]);
                 
  void
  MACGradient(BoxData<double, DIM>   a_macGrad[DIM],
              BoxData<double, 1  >&  a_phicc);
                 
  void
  averageGradientToCell(BoxData<double, DIM>& a_gradcell,
                        BoxData<double, DIM>  a_gradface[DIM]);
  void
  solveElliptic(BoxData<double, 1> & a_phi,
                BoxData<double, 1> & a_rhs,
                const double& a_alpha, const double& a_beta, 
                const string& a_basename);


  void
  getNuLaplU(BoxData<double, 1> & a_source, 
             BoxData<double, 1> & a_scalarCell,
             double a_coeff);

  Box    m_domain;
  double m_dx;
  double m_viscosity;
  int    m_nghost;

  void  defineStencils();

  Stencil<double> m_velToFaceSten[DIM];
  Stencil<double> m_gradToCellSten[DIM];
  Stencil<double> m_macDivergeSten[DIM];
  Stencil<double> m_macGradientSten[DIM][DIM];

#if DIM==2
  Stencil<double> m_vortsten[2];
#else
  Stencil<double> m_vortsten[3][3];
#endif

};
#endif //end include guard
