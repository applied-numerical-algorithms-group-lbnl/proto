PROTO_KERNEL_START
void
RKTaylorCoefs(Var<double,,RKORDER*NUMCOMPS>& a_UTaylor,
              const Var<double,NUMCOMPS>& a_kstage,
              const Var<double,NUMCOMPS>& a_U,
              int a_stage
              )
{
  for (int n = 0; n < NUMCOMPS;n++)
    {
      if (a_stage == 0) 
        {
          a_UTaylor(n*NUMCOMPS) = a_U(n);
        }
      if (a_stage == 1)
        {
          a_UTaylor(1 + n*NUMCOMPS) = a_kstage(n);
          a_UTaylor(2 + n*NUMCOMPS) = -(3.0/2)*a_kstage(n);
          a_UTaylor(3 + n*NUMCOMPS) = (2.0/3)*a_kstage(n);
        }
      if (a_stage == 2)
        {
          a_UTaylor(2 + n*NUMCOMPS) += a_kstage(n);
          a_UTaylor(3 + n*NUMCOMPS) += (-2.0/3)*a_kstage(n);
          
        }
      if (a_stage == 3)
        {
          a_UTaylor(2 + n*NUMCOMPS) += (-1.0/2)*a_kstage(n);
          a_UTaylor(3 + n*NUMCOMPS) += (2.0/3)*a_kstage(n);
        }
    }
}
PROTO_KERNEL_END(RKTaylorCoefsF, RKTaylorCoefs)
PROTO_KERNEL_START
void
RKTIimeInterp(Var<double,NUMCOMPS>& a_U,
              Var<double,RKORDER*NUMCOMPS>& a_UTaylor,
              double& a_chi
              )
{
  for (int n = 0; n < NUMCOMPS;n++)
    {
      a_U(n) = 0.;
      for (int term = 0; term < RKORDER; term++)
        {
          a_U(n) = a_U(n)*a_chi + a_UTaylor(RKORDER-term-1 + n*NUMCOMPS);
        }
    }
}
PROTO_KERNEL_END(RKTimeInterpF, RKTimeInterp)
template
<class OpType, typename T,unsigned int C=1,Memtype MEM>
void AMRRKExplicit<OpType,T,C,MEM>::define
                                   (shared_ptr<AMRData<T,C,MEM> > a_dataPtr,
                                    const double& a_dx,
                                    int a_timeRefRatio,
                                    int a_level,
                                    )
{
  // Set up simple data members.
  m_dataPtr = a_dataPtr;
  m_level = a_level;
  m_dbl = m_dataPtr.grids()[m_level];
  m_refRatio = PR_AMR_REFRATIO*Point::Ones();
  m_timeRefRatio = a_timeRefRatio;
  m_dx = a_dx;

  // define m_op, get information about the size of its ghost cell region.
  m_op.define(m_dx);
  m_ghost = m_op.ghostsize();
  
  // If we have a finer level, make recursive call to constructor, set up storage
  // for the saving the Taylor expansion coefficients for the current level. 
  if (m_level < m_dataPtr.size() - 1)
    {
      dxFiner = m_dx/PR_AMR_REFRATIO;
      m_finer = shared_ptr<AMRRKExplicit<OpType,T,C,MEM> >
                            (new AMRRKExplicit<OpType,T,C,MEM>
                             (m_dataPtr,
                              dxFiner,
                              m_timeRefRatio,
                              m_level + 1)
                             );
     
      m_UTaylor.define(a_dataPtr.grid()[m_level],Point::Zeros());
    }
  // If we have a coarser level:
  // (1) build interpolation stencil,
  // (2) set up storage for Taylor expansion coefficients from the next coarser level,
  // (3) For each patch, compute boxes over which the ghost cell interpolation will
  // be applied.
  if (m_level > 0)
    {
      // Set up for fourth-order in space and time c/f bcs.

      // Storage for RK4 interpolation function
      Point ghostInterp = m_op.ghostSize()/m_refRatio + (3*Point::Ones());
      M_UTaylorCoarsened.define(
                                a_dataPtr.grid()[m_level].coarsen(m_refRatio),
                                ghostInterp);
      // Construct fourth-order spatial stencil.
      Box bStencil(*Point::Ones(-1),Point::Ones());
      vector<Point> stencilPts;
                                                         
      for (auto bit = bStencil.begin();!bit.end();++bit)
        {
          if (*bit!=Point::Zeros()) stencilPts.push_back(*bit);
        }
      for (SideIterator sit.begin();!sit.done();++sit)
        {
          for (int dir = 0;dir < DIM; dir++)
            {
              stencilPts.push_back(Point::Basis(dir,sit().sign()));              
            }
        }
      m_cfInterp.Build(stencilPts,4,PR_AMR_RATIO);
      
      // Vector of boxes to which to interpolate c/f bcs.
      
      m_adjacentBoxes.resize(m_dbl.mySize());
      for (dit = m_dbl.begin(); *dit!=dit.end();++dit)
        {
          int myindex = m_dbl.myIntIndex(*dit);
          Box myBox = m_dbl[*dit].grow(m_dataPtr[m_level].ghost());          
          vector<Box>& myAdjBoxes = m_adjacentBoxes[myIndex];
          Box bx(Point::Ones(-1),Point::Ones());
          
          // We are assuming here that the ghost cell region is no wider than a patch,
          // so that the ghost cells on the same level can be filled from patches in
          // the 3^DIM cube of adjacent patches.
          
          Point ptMyBox = m_dbl.point(*dit);
          
          for (auto bxit = bx.begin();*bit!=bx.end();++bxit)
            {
              if (*bxit != Point::Zeros())
                {
                  Point pt = ptMyBox + *bxit;
                  Box bxPt = Box(pt,pt).refine(m_dbl.boxsize());
                  if (m_dbl.find(pt) == m_dbl.end())
                    {
                      myAdjBoxes.push_back(myBox&bxPt);
                    }
                }
            }
        }     
    }
}
template
<class OpType, typename T,unsigned int C=1,Memtype MEM>
void AMRRKExplicit<OpType,T,C,MEM>::advance
                                   (LevelFluxRegister<T,C,MEM> a_coarseLFR,
                                    const double& a_dt,
                                    const double& a_chi,
                                    bool a_willRegridCoarse)
{
  // Setup for interpolation from coarse grid data using dense representation.
  
  LevelBoxData<T,C,MEM> UStage(m_amrgrids[m_level],m_ghost);
  (*m_dataPtr)[m_level].copyTo(UStage);
  
  for (int stage = 0; stage < RKSTAGES; stage++) // Stage loop.
    {
      if (stage != 0) UStage.exchange();
      
      for (auto dit = m_dbl.begin();dit != dit.end();++dit) 
        {
          // interpolate coarse grid data in time, then in space, to fill ghost cells.
          auto & ULocal = UStage[*dit];
          Array<BoxData<T,C,MEM>, DIM> fluxes;
          
          for (int dir = 0; dir < DIM; dir++)
            {
              fluxes[dir].define(XXX);
            }
          chiStage = a_chi + stage*m_stageTimes/m_timeRefRatio;
          interpCoarse(ULocal,UTaylorLocal[*dit],chiStage);  
              
          // Perform stage calculation.

          m_op(ULocal,kStage,fluxes,a_dt*m_stageWgt[stage]);             
              
          if (m_level < m_dataPtr.size()-1)
            {
              // Store stage information.
              forallInPlace(RKTaylorCoefs,m_UTaylor[*dit],kStage,stage);
                  
              // Increment flux registers.                
              for (int dir = 0; dir < DIM; dir++)
                {
                  m_register.incrementCoarse(fluxes[dir],*dit,-a_dt*m_wgt[stage],dir);
                }
            }
          if (m_level> 0)
            {
              // Increment fine flux registers.
              for (int dir = 0; dir < DIM; dir++)
                {
                  a_coarseRegister.incrementFine(fluxes[dir],*dit,a_dt*m_wgt[stage],dir);
                }
            }
        } // end loop for a single stage.
    } // end loop over RK4 stages.

  // recursive call to advance finer levels.

  // Copy Taylor coefs to the data holder on the next finer level.
  m_UTaylor.copyTo(m_finer->m_UTaylorCoarse);

  // Advance finer levels.
  m_numSteps++;
  bool willRegrid = (m_numSteps == m_regridInterval - 1);
  for (int step = 0; step < m_timeRefRatio;step++)
    {
      double chi = step*1.0/m_timeRefRatio;
      double dtFine = dt/m_timeRefRatio;
      m_finer->advanceLevel(m_register,
                            dtFine,
                            chi,
                            willRegrid);
    }
  
  // Reflux.
  double dxI = 1.0/m_dx;
  m_register.reflux((*m_dataPtr)[m_level],dxI);
  
  // Regrid. Regridding will be managed by the finest level that doesn't change.
  if (m_level < (*dataPtr).size() - 1)&&(m_numSteps == m_regridInterval-1))
    {
      m_numSteps = 0;
      if (!a_willRegridCoarse)
        {
          // call regridding - TBD
          cout << "Warning - you've hit regridding which is not yet implemented" << endl;
        }
    }
}
