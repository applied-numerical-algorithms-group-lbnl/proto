template<typename T, MemType MEM>
Multigrid<T,MEM>::Multigrid(
                            DisjointBoxLayout a_finestLayout,
                            int a_solverType,
                            int a_numLevels)
{
  m_p = 2*DIM;
  m_pbottom = 10;
  m_numLevels = a_numLevels;
  m_smoother = Stencil<T>::Laplacian();

  m_laplacian = LaplaceStencil<T>(a_solverType);
  if (procID() == 0) m_laplacian.print();
  m_order = 2;
  if (a_solverType > 0) m_order = 4;
  m_ghost = m_laplacian.ghost()[0];
  m_layouts.resize(a_numLevels+1);
  if (procID() == 0)
    {
      cout << "input layout " << endl;
      a_finestLayout.print();
      cout << a_solverType << ", " << a_numLevels << endl;
    }
  Point fineboxsize = a_finestLayout.boxSize();
  DisjointBoxLayout curlayout = a_finestLayout;
  // construct dbls for all levels.
  // cout << procID() << ", " << a_numLevels << endl;

  for (int lev = a_numLevels; lev > -1 ; lev--)
    {    
      m_layouts[lev] = curlayout;
      if (lev == 0) break;
      ProblemDomain curdomain = curlayout.domain();
      ProblemDomain coarsedomain = curdomain.coarsen(Point::Ones(2));
      if (coarsedomain.coarsenable(fineboxsize))
        {
          curlayout = DisjointBoxLayout(coarsedomain, fineboxsize);
        }
      else
        {
          auto coarsesize =  coarsedomain.box().high() + Point::Ones();
          curlayout = DisjointBoxLayout(coarsedomain,coarsesize);
        }
    }

  // Extract diagonal coefficient for the stencils.
  for (int index = 0; index < m_smoother.size(); index++) 
    {
      if (m_smoother.offsets()[index] == Point::Zeros())
        {
          m_diagCoefSmooth = m_smoother.coefs()[index];
        }
    }
    for (int index = 0; index < m_laplacian.size(); index++) 
    {
      if (m_laplacian.offsets()[index] == Point::Zeros())
        {
          m_diagCoefLaplacian =  m_laplacian.coefs()[index];
        }
    }
    if (procID() == 0) cout
                         << "diag coef laplacian = "
                         <<  m_diagCoefLaplacian << endl;
    
    //#if 0
    for (int dir = 0; dir < DIM; dir++)
      {
        Point bv = Point::Basis(dir,-1);
        m_sten2low[dir] = (-1.0)*Shift(-bv);       
        m_sten42low[dir] = (-18.0)*Shift(-bv) +8.0*Shift(-2*bv) +
          (-9.0/5.0)*Shift(-3*bv);
        bv = Point::Basis(dir,1);
        m_sten2high[dir] = (-1.0)*Shift(-bv);      
        m_sten42high[dir] = (-18.0)*Shift(-bv) +8.0*Shift(-2*bv) +
          (-9.0/5.0)*Shift(-3*bv);
      }

    Box kernel = Box::Kernel(1);
    m_sten41ij.define(kernel);
    int kk = 0;
    m_sten4.clear();
    for (auto bit : kernel)
      {
        if (bit != Point::Zeros())
          {
            
            Point bit2 = bit;
            Stencil<T> sten41bit =
              (    -3.0)*Shift(-bit2) +
              (     1.0)*Shift(-2*bit2) +
              (-1.0/5.0)*Shift(-3*bit2);
            //Stencil<T> sten41bit = -1.*Shift(-bit2);
            m_sten4.push_back(sten41bit);
            m_sten41ij(bit2) = kk;
            kk++;
          }
      }

    //#endif
}
template<typename T, MemType MEM>
void Multigrid<T,MEM>::mgRelax(LevelBoxData<T,1,MEM>& a_phi,
                               const LevelBoxData<T,1,MEM>& a_rhs,
                               T a_h,
                               int a_level)
{
  PR_TIMERS("Multigrid::mgRelax");
  HDF5Handler h5;
  if (a_level > 0)
    {
      pointRelax(a_phi,a_rhs,a_h,m_p,a_level);
      
      auto layout = m_layouts[a_level];
      auto layoutC = m_layouts[a_level-1];
      auto layoutCoarsened = layout.coarsen(Point::Ones(2));
      LevelBoxData<T,1,MEM> phiCoarse(layoutC,Point::Ones());      
      LevelBoxData<T,1,MEM> phiCoarsened(layoutCoarsened,Point::Ones());   
      LevelBoxData<T,1,MEM> resCoarsened(layoutCoarsened,Point::Zeros());
      LevelBoxData<T,1,MEM> phiCoarsened0(layoutCoarsened,Point::Ones());
      LevelBoxData<T,1,MEM> rhsCoarse(layoutC,Point::Zeros());
      auto pd = layout.domain();
      Stencil<T> op = m_smoother;     
      int ghost = op.ghost()[0];
      a_phi.exchange();
      for (auto dit : layout)
        {
          BoxData<T,1,MEM> res(layout[dit]);
          setGhost(a_phi[dit],a_h,pd,layout[dit]);
          res |= op(a_phi[dit],-1.0/(a_h*a_h));  
          res += a_rhs[dit];
          resCoarsened[dit] |= Stencil<T>::AvgDown(2)(res);
          phiCoarsened[dit] |= Stencil<T>::AvgDown(2)(a_phi[dit]);
        }
      phiCoarsened.copyTo(phiCoarse);
      phiCoarsened.copyTo(phiCoarsened0);
      resCoarsened.copyTo(rhsCoarse);
      auto pdC = layoutC.domain();
      phiCoarse.exchange();
      for (auto dit : layoutC)
        {
          setGhost(phiCoarse[dit],2*a_h,pdC,layoutC[dit]);
          rhsCoarse[dit] += m_smoother(phiCoarse[dit],1.0/(4*a_h*a_h));
        }
      mgRelax(phiCoarse,rhsCoarse,2*a_h,a_level-1);      
      phiCoarse.copyTo(phiCoarsened);
      
      for (auto dit : layout)
        {
          phiCoarsened[dit] -= phiCoarsened0[dit];
          a_phi[dit] +=
            InterpStencil<T>::Constant(Point::Ones(2))(phiCoarsened[dit]);
        }    
      pointRelax(a_phi,a_rhs,a_h,m_p,a_level);
    }
  else
    {
      pointRelax(a_phi,a_rhs,a_h,m_pbottom,a_level);
    }
}
template<typename T, MemType MEM>
void
Multigrid<T,MEM>::residual(
                           LevelBoxData<T,1,MEM>& a_residual,
                           LevelBoxData<T,1,MEM>& a_phi,
                           const LevelBoxData<T,1,MEM>& a_rhs,
                           T a_h)
{
  // compute residual of opertor you are solving for.
  PR_TIMERS("Multigrid<T,MEM>::residual");
  auto layout = a_phi.layout();
  auto pd = layout.domain();
  a_phi.exchange();
  int ghost = m_laplacian.ghost()[0];
  for (auto dit : layout)
    {
      setGhost(a_phi[dit],a_h,pd,layout[dit],m_order,ghost);
      a_residual[dit] |= m_laplacian(a_phi[dit],-1./(a_h*a_h));
      a_residual[dit] += a_rhs[dit];
    }
}
template<typename T, MemType MEM>
void
Multigrid<T,MEM>::pointRelax(LevelBoxData<T,1,MEM>& a_phi, 
                             const LevelBoxData<T,1,MEM>& a_rhs,
                             T a_h, 
                             int a_numIter,
                             int a_level)
{
  PR_TIMERS("Multigrid<T,MEM>::pointRelax");
  T rcoef = -a_h*a_h/(m_diagCoefSmooth*2);
  auto layout = m_layouts[a_level];
  auto pd = a_rhs.layout().domain();
  for (int iter =  0; iter < a_numIter; iter++)
    {
      a_phi.exchange();
      for (auto dit : layout)
        {
          setGhost(a_phi[dit],a_h,pd,layout[dit]);
          BoxData<T,1,MEM> LOfPhi = m_smoother(a_phi[dit],1.0/(a_h*a_h));  
          forallInPlace([] PROTO_LAMBDA(Var<T>& a_phi,
                                        const Var<T>& a_LOfPhi,
                                        const Var<T>& a_rhs,
                                        T a_rcoef)
                        {
                          a_phi(0) += (a_LOfPhi(0) - a_rhs(0))*a_rcoef;
                        },a_phi[dit],LOfPhi,a_rhs[dit],rcoef);
        }
    }
}
template<typename T, MemType MEM>
void Multigrid<T,MEM>::setGhost(BoxData<T,1,MEM>& a_phi,
                                T a_h,
                                ProblemDomain a_domain,
                                Box a_bx,
                                int a_order,
                                int a_ghost
                                )
{
  PR_TIMERS("Multigrid<T,MEM>::setGhost");
  HDF5Handler h5;
  Box domainbox = a_domain.box();
  //if (a_order == 2)
    {
      PROTO_ASSERT(a_ghost == 1,"setGhost assert failure");
      for (int dir = 0; dir < DIM; dir++)
        {
          Point bv = Point::Basis(dir,-1);
          Box bxlo = a_bx.adjacent(bv);
          if (!domainbox.contains(bxlo))
            {
              BoxData<T,1,MEM> phiShift(bxlo);
              phiShift |= m_sten2low[dir](a_phi);
              forallInPlace_p(f_bc,bxlo,a_phi,phiShift,a_h,m_dirichlet2,bv);
            }
          bv = Point::Basis(dir,1);
          Box bxhi =a_bx.adjacent(bv);
          if (!domainbox.contains(bxhi))
            {
              BoxData<T,1,MEM> phiShift(bxhi);
              phiShift |= m_sten2high[dir](a_phi);
              forallInPlace_p(f_bc,bxhi,a_phi,phiShift,a_h,m_dirichlet2,bv);
            }
        }
    }
  return;
    if (a_order == 4)
    {
      for (auto bit : Box::Kernel(1))
        {
          Point bit2 = bit;
          if (bit2 != Point::Zeros())
            {
              int index = m_sten41ij(bit2);
              Box bxshift = a_bx.adjacent(bit2);
              if (!domainbox.contains(bxshift))
                {
                  BoxData<T,1,MEM> phiShift(bxshift);
                  phiShift |= m_sten4[index](a_phi); 
                  forallInPlace_p
                    (f_bc,bxshift,a_phi,phiShift,a_h,m_dirichlet41,bit);
                }
            }
        }
      if (a_ghost == 2)
        {
          abort();
          for (int dir = 0; dir < DIM; dir++)
            {
              Point bv = Point::Basis(dir,-2);
              Box bxlo = a_bx.adjacent(bv);
              if (!domainbox.contains(bxlo))
                {
                  bxlo.high()[dir] -= 1;
                  BoxData<T,1,MEM> phiShift(bxlo);
                  phiShift |= m_sten42low[dir](a_phi);
                  forallInPlace_p
                    (f_bc,bxlo,a_phi,phiShift,a_h,m_dirichlet42,bv);
                }
              bv = Point::Basis(dir,2);
              Box bxhi =a_bx.adjacent(bv); 
              if (!domainbox.contains(bxhi))
                {
                  bxhi.low()[dir] += 1;
                  BoxData<T,1,MEM> phiShift(bxhi);
                  phiShift |= m_sten42high[dir](a_phi);
                  forallInPlace_p
                    (f_bc,bxhi,a_phi,phiShift,a_h,m_dirichlet42,bv);
                }
            }
        }
    }

}
template<typename T, MemType MEM>
T Multigrid<T,MEM>::resnorm(LevelBoxData<T,1,MEM>& a_phi,
                            const LevelBoxData<T,1,MEM>& a_rho,
                            T a_h)
{
  Reduction<T,Abs,MEM> rxn;
  rxn.reset();
  auto layout = a_phi.layout();
  auto pd = layout.domain();
  auto ghost = m_laplacian.ghost()[0];
  for (auto dit : layout)
    {
      setGhost(a_phi[dit],a_h,pd,layout[dit],m_order,ghost);
      BoxData<T,1,MEM> res(layout[dit]);
      res |= m_laplacian(a_phi[dit],-1.0/(a_h*a_h));
      res += a_rho[dit];
      T maxres = res.absMax();
      rxn.reduce(&maxres,1);
    }  
  T retval = rxn.fetch();
  return retval;
}
