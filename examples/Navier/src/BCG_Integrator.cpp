#include "BCG_Integrator.H"
#include "SGMultigrid.H"
#include "GodunovAdvectionOp.H"
#include "CommonTemplates.H"
#include "Proto.H"
#include "Proto_Timer.H"
using std::cout;
using namespace Proto;
void
bcgEnforceVecBC(BoxData<double, DIM>& a_phi, 
                int a_numghost)
{
  for(int idir = 0; idir < DIM; idir++)
  {
    protocommon::enforceSGBoundaryConditions<double, DIM>(a_phi, a_numghost, idir);
  }
};
void
bcgEnforceScaBC(BoxData<double, 1>& a_phi, 
                int a_numghost)
{
  for(int idir = 0; idir < DIM; idir++)
  {
    protocommon::enforceSGBoundaryConditions<double, 1>(a_phi, a_numghost, idir);
  }
};
////
void
BCG_Integrator::
averageVelocityToFaces(BoxData<double, 1  >   a_velface[DIM],
                       BoxData<double, DIM> & a_velcell)
{
  PR_TIME("bcg:aveveltofaces");
  bcgEnforceVecBC(a_velcell, m_nghost);
  for(int idir = 0; idir < DIM; idir++)
  {
    const Box& facebx =  m_domain.extrude(idir);
    BoxData<double, 1> velcomp  = slice(a_velcell, idir);
    a_velface[idir] |= m_velToFaceSten[idir](velcomp, facebx);
  }
}
////
void
BCG_Integrator::
averageGradientToCell(BoxData<double, DIM>   & a_gradcell,
                      BoxData<double, DIM>     a_gradface[DIM])
{
  PR_TIME("bcg:avegradtocells");
  for(int idir = 0; idir < DIM; idir++)
  {
    BoxData<double, 1> facecomp  = slice(a_gradface[idir], idir);
    BoxData<double, 1> cellcomp  = slice(a_gradcell      , idir);
    cellcomp  |= m_gradToCellSten[idir](facecomp, m_domain);
  }
}
////
void
BCG_Integrator::
MACDivergence(BoxData<double, 1> & a_divergence,
              BoxData<double, 1>   a_velface[DIM])
{
  PR_TIME("bcg:macdiv");
  a_divergence.setVal(0.);
  for(int idir = 0; idir < DIM; idir++)
  {
    a_divergence += m_macDivergeSten[idir](a_velface[idir], m_domain);
  }
}


////
void
BCG_Integrator::
MACGradient(BoxData<double, DIM>   a_macGrad[DIM],
            BoxData<double, 1  >&  a_phicc)
{
  PR_TIME("bcg:macgrad");
  for(int faceDir = 0; faceDir < DIM; faceDir++)
  {
    for(int gradDir = 0; gradDir < DIM; gradDir++)
    {
      const Box& facebx =  m_domain.extrude(faceDir);
      BoxData<double, 1> gradcomp  = slice(a_macGrad[faceDir], gradDir);
      gradcomp |= m_macGradientSten[faceDir][gradDir](a_phicc, facebx);
    }
  }
}

//divides face centered velocity into vortical compoenent and a gradient of a scalar
void
BCG_Integrator::
MACProject(BoxData<double,   1> a_velocity[DIM],
           BoxData<double, DIM> a_gradpres[DIM])
{
  PR_TIME("bcg:macproject");
  Box grownBox = m_domain.grow(1);
  BoxData<double, 1> divergence(m_domain);
  BoxData<double, 1>     scalar(grownBox);
  MACDivergence(divergence, a_velocity);
  double alpha = 0; double beta = 1; //solving Poisson
  double sumdiv = divergence.sum();
  std::cout << "sum of divergence = " << sumdiv << std::endl;
  sumdiv /= m_domain.size();
  divergence -= sumdiv;
  solveElliptic(scalar, divergence, alpha, beta, string("projection:"));
//  BoxData<double, DIM> gradient[DIM];
//  for(int idir = 0; idir < DIM; idir++)
//  {
//    const Box& facebx =  m_domain.extrude(idir);
//    gradient[idir].define(facebx);
//  }
  MACGradient(a_gradpres, scalar);
  for(int idir = 0; idir < DIM; idir++)
  {
    BoxData<double, 1> gradcomp  = slice(a_gradpres[idir], idir);
    a_velocity[idir] -= gradcomp;
  }  
}
//divides cell centered velocity into vortical compoenent and a gradient of a scalar
void
BCG_Integrator::
ccProject(BoxData<double, DIM>& a_velocity,
          BoxData<double, DIM>& a_gradpres)
{
  PR_TIME("bcg:ccproject");
  BoxData<double, 1  > faceVel[DIM];
  BoxData<double, DIM> faceGrad[DIM];
  for(int idir = 0; idir < DIM; idir++)
  {
    const Box& facebx =  m_domain.extrude(idir);
    faceVel[idir ].define(facebx);
    faceGrad[idir].define(facebx);
  }
  averageVelocityToFaces(faceVel, a_velocity);
  MACProject(faceVel, faceGrad);
  averageGradientToCell(a_gradpres, faceGrad);
  a_velocity -= a_gradpres;
}
////
void
BCG_Integrator::
getNuLaplU(BoxData<double, 1> & a_source, 
           BoxData<double, 1> & a_scalarCell,
           double a_coeff)
{
  PR_TIME("bcg::nulapphi");
  //compue nu*lapl u
  double alpha = 0; double beta = a_coeff;
  SGMultigrid solver(alpha, beta, m_dx, m_domain);
  solver.applyOp(a_source, a_scalarCell);
}
////
void
BCG_Integrator::
getUDotDelU(BoxData<double, DIM> & a_udelu,
            BoxData<double, DIM> & a_velocity,
            const double         & a_dt)
{
  PR_TIME("bcg::getudelu");
  //using velocity averaged to faces as advection velocity
  //just used for upwinding
  BoxData<double, 1>   advectVel[DIM];
  BoxData<double, DIM> faceGrad[DIM];
  BoxData<double, DIM> faceVelo[DIM];
  for(int idir = 0; idir < DIM; idir++)
  {
    const Box& facebx =  m_domain.extrude(idir);
    faceGrad[idir ].define(facebx);
    advectVel[idir].define(facebx);
    faceVelo[idir ].define(facebx);
  }
  averageVelocityToFaces(advectVel, a_velocity);

  int ideb = 0;
  GodunovAdvectionOp advOp(m_dx);

  int doingvel = 1;
  for(int velComp = 0; velComp < DIM; velComp++)
  {
    BoxData<double,1>  scalarFace[DIM];
    BoxData<double,1>  scalarCell = slice(a_velocity, velComp);
    for(int faceDir = 0; faceDir < DIM; faceDir++)
    {
      scalarFace[faceDir] = slice(faceVelo[faceDir], velComp);
    }
    BoxData<double, 1> source(m_domain);
    
    getNuLaplU(source, scalarCell, m_viscosity);
    bcgEnforceScaBC(scalarCell, m_nghost);
    advOp.advectToFaces(scalarFace, advectVel, scalarCell, 
                        source, a_velocity, m_domain, doingvel, a_dt);
    ideb++;
    
  }

  //now need to mac project -- reuse advVel to get the proper pressure gradient
  //because we need a scalar velocity
  for(int faceDir = 0; faceDir < DIM; faceDir++)
  {
    BoxData<double, 1> velcomp = slice(faceVelo[faceDir], faceDir);
    velcomp.copyTo(advectVel[faceDir]);
  }
  MACProject(advectVel, faceGrad);
  //now subtract the gradient of the scalar off all components 
  //(advectvel was corrected in mac project and just contains normal compoenents)
  //we use advectvel again below
  for(int faceDir = 0; faceDir < DIM; faceDir++)
  {
    faceVelo[faceDir] -= faceGrad[faceDir];
  }

  //compute the divergence of the flux
  for(int velComp = 0; velComp < DIM; velComp++)
  {
    BoxData<double,1>  velcomp[DIM];
    for(int faceDir = 0; faceDir < DIM; faceDir++)
    {
      velcomp[faceDir] = slice(faceVelo[faceDir], velComp);
    }
    BoxData<double, 1> divFcomp = slice(a_udelu, velComp);
    advOp.divUPhi(divFcomp, advectVel, velcomp, m_domain);
  }
}
////
void
BCG_Integrator::
advanceSolution(BoxData<double, DIM>& a_velocity, 
                BoxData<double, DIM>& a_gradpres, 
                const double        & a_dt)
{
  PR_TIME("bcg::advancesolution");
  //get advective derivative
  BoxData<double, DIM> udelu(m_domain);
  getUDotDelU(udelu, a_velocity, a_dt);

  BoxData<double, DIM> ustar(a_velocity.box());
  if(m_viscosity > 1.0e-16)
  {
    //get 0.5*nu*lapl(vel) for crank-nicolson rhs
    BoxData<double, DIM> halfnulaplu(m_domain);
    for(int velComp = 0; velComp < DIM; velComp++)
    {
      BoxData<double,1>  velcomp = slice(a_velocity, velComp);
      BoxData<double,1>  lapcomp = slice(halfnulaplu, velComp);
      getNuLaplU(lapcomp, velcomp, 0.5*m_viscosity);
    }

    //form c-n rhs
    BoxData<double, DIM> rhs(m_domain);
    rhs.setVal(0.);
    rhs -= udelu;
    rhs -= a_gradpres;
    rhs += halfnulaplu;
    rhs *= a_dt;
    rhs += a_velocity;

    //solve for u*
    double alpha = 1;
    double beta =  -0.5*m_viscosity*a_dt;
    for(int velComp = 0; velComp < DIM; velComp++)
    {
      BoxData<double,1>  velcomp = slice(ustar, velComp);
      BoxData<double,1>  rhscomp = slice(rhs,   velComp);
      solveElliptic(velcomp, rhscomp, alpha, beta, string("viscous_solve:"));
    }

  }
  else
  {
    ustar.setVal(0.);
    ustar -= udelu;
    ustar += a_gradpres;
//    ustar -= a_gradpres;
    ustar *= a_dt;
    ustar += a_velocity;
  }
  //now we need to project the solution onto its divergence-free subspace.
  //we include the usual tricks to improve stability
  //w = vel + dt*gph^n-1/2
  //u = P(w)
  //grad p^n+1/2 = (1/dt)(I-P)w;
  ccProject(ustar, a_gradpres);
  ustar.copyTo(a_velocity);
  a_gradpres /= a_dt;

}
///
void
BCG_Integrator::
solveElliptic(BoxData<double, 1> & a_phi,
              BoxData<double, 1> & a_rhs,
              const   double     & a_alpha, 
              const   double     & a_beta,
              const string& a_solvename)
{
  PR_TIME("bcg::solveelliptic");
  using std::cout;
  using std::endl;
  int numsmooth = 4;
  int maxiter = 27;
  double tol = 1.0e-9;
  SGMultigrid solver(a_alpha, a_beta, m_dx, m_domain);
  SGMultigrid::s_numSmoothUp   = numsmooth;
  SGMultigrid::s_numSmoothDown = numsmooth;

  BoxData<double, 1> res(m_domain);

  a_phi.setVal(0.);
  int iter = 0;
  double rhsmax = a_rhs.absMax();
  double resStart = std::max(a_rhs.absMax(), tol);
  double resIter  = rhsmax;
  double resIterOld = resIter;
  cout << a_solvename << "iter = " << iter << ", ||resid|| = " << resIter << endl;
  double mintol = 1.0e-10;
  double ratio = 2;
  double minrat = 1.05;
  while((resIter > mintol) && (resIter > tol*resStart) && (iter <  maxiter) && (ratio >minrat))
  {
    solver.vCycle(a_phi, a_rhs);
    solver.residual(res, a_phi, a_rhs);
  
    iter++;
    resIter = res.absMax();
    ratio = resIterOld/resIter;
    resIterOld = resIter;
    cout << "iter = " << iter << ", ||resid|| = " << resIter << endl;
  }
  solver.enforceBoundaryConditions(a_phi);
}
////
void
BCG_Integrator::
defineStencils()
{
  PR_TIME("bcg::definestencils");
  double dx = m_dx;
  //stencils to average velocities from cells to faces and to increment the divergence 
  // of a mac velocity
  for(int faceDir = 0; faceDir < DIM; faceDir++)
  {
    m_velToFaceSten[faceDir]  = (0.5)*Shift::Zeros() +  (0.5)*Shift::Basis(faceDir,-1.0);    
    m_gradToCellSten[faceDir] = (0.5)*Shift::Zeros() +  (0.5)*Shift::Basis(faceDir, 1.0);    
    m_macDivergeSten[faceDir] = (-1.0/dx)*Shift::Zeros() + (1.0/dx)*Shift::Basis(faceDir, 1.0);    
    m_macGradientSten[faceDir][faceDir] = (1.0/dx)*Shift::Zeros() + (-1.0/dx)*Shift::Basis(faceDir, -1.0);    
    for(int gradDir = 0; gradDir < DIM; gradDir++)
    {
      if(faceDir != gradDir)
      {
        Point hihi = Point::Basis(gradDir, 1.0);
        Point hilo = Point::Basis(gradDir,-1.0);
        Point lolo = Point::Basis(gradDir,-1.0) + Point::Basis(faceDir, -1);
        Point lohi = Point::Basis(gradDir, 1.0) + Point::Basis(faceDir, -1);
        m_macGradientSten[faceDir][gradDir] = 
          ( 0.25/dx)*Shift(hihi) + ( 0.25/dx)*Shift(hilo)  +
          (-0.25/dx)*Shift(lolo) + (-0.25/dx)*Shift(lohi);
      }
    }
  }

///vorticity stencils
#if DIM==2
  m_vortsten[0] = (-0.5/dx)*Shift::Basis(1, 1.0) +  (0.5/dx)*Shift::Basis(1, -1.0);
  m_vortsten[1] =  (0.5/dx)*Shift::Basis(0, 1.0) + (-0.5/dx)*Shift::Basis(0, -1.0);
#else
  m_vortsten[0][2] = (-0.5/dx)*Shift::Basis(1, 1.0) +  (0.5/dx)*Shift::Basis(1, -1.0);
  m_vortsten[1][2] =  (0.5/dx)*Shift::Basis(0, 1.0) + (-0.5/dx)*Shift::Basis(0, -1.0);
  m_vortsten[2][2] = (0.0)*Shift::Zeros();

  m_vortsten[0][1] = (-0.5/dx)*Shift::Basis(2, 1.0) +  (0.5/dx)*Shift::Basis(2, -1.0);
  m_vortsten[2][1] =  (0.5/dx)*Shift::Basis(0, 1.0) + (-0.5/dx)*Shift::Basis(0, -1.0);
  m_vortsten[1][1] = (0.0)*Shift::Zeros();

  m_vortsten[2][0] =  (0.5/dx)*Shift::Basis(1, 1.0) + (-0.5/dx)*Shift::Basis(1, -1.0);
  m_vortsten[1][0] =  (0.5/dx)*Shift::Basis(2, 1.0) +  (0.5/dx)*Shift::Basis(2, -1.0);
  m_vortsten[0][0] = (0.0)*Shift::Zeros();
#endif
}

//cheerfully stolen from the euler example
void
BCG_Integrator::
enforceBoundaryConditions(BoxData<double, DIM>& a_vel)
{
  PR_TIME("bcg::enforcebcs");
  for(int idir = 0; idir < DIM; idir++)
  {
    Point dirlo = -Point::Basis(idir);
    Point dirhi =  Point::Basis(idir);
    Box dstbxlo = m_domain.edge(dirlo);
    Box dstbxhi = m_domain.edge(dirhi);

    //these are swapped because you have to move in different 
    //directions to get the periodic image
    Point shifthi = dirlo*m_domain.size(idir);
    Point shiftlo = dirhi*m_domain.size(idir);

    Box srcbxlo = dstbxlo.shift(shiftlo);
    Box srcbxhi = dstbxhi.shift(shifthi);

    a_vel.copy(a_vel, srcbxlo, 0, dstbxlo, 0, DIM);
    a_vel.copy(a_vel, srcbxhi, 0, dstbxhi, 0, DIM);
  }
}
#if DIM==2
///
void 
BCG_Integrator::
computeVorticity(BoxData<double, 1>         & a_vorticity,
                 BoxData<double, DIM>       & a_velocity)
{
  PR_TIME("bcg::vorticity");
  enforceBoundaryConditions(a_velocity);
  a_vorticity.setVal(0.);
  for(int idir = 0; idir < DIM; idir++)
  {
    BoxData<double, 1> velcomp  = slice(a_velocity, idir);
    a_vorticity += m_vortsten[idir](velcomp, m_domain);
  }
}
//
#else
void 
BCG_Integrator::
computeVorticity(BoxData<double, DIM>       & a_vorticity,
                 BoxData<double, DIM>       & a_velocity)
{
  PR_TIME("bcg::vorticity3d");
  enforceBoundaryConditions(a_velocity);
  a_vorticity.setVal(0.);
  for(int veldir = 0; veldir < DIM; veldir++)
  {
    for(int vortdir = 0; vortdir < DIM; vortdir++)
    {
      if(veldir != vortdir)
      {
        BoxData<double, 1> velcomp   = slice(a_velocity , veldir);
        BoxData<double, 1> vortcomp  = slice(a_vorticity, vortdir);
        vortcomp += m_vortsten[veldir][vortdir](velcomp, m_domain);
      }
    }
  }
}
#endif
