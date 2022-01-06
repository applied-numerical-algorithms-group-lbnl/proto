#include "AdvectionTest.H"
PROTO_KERNEL_START
void
RKTaylorCoefsF(Var<double,RKORDER*NUMCOMPS>& a_UTaylor,
              Var<double,NUMCOMPS>& a_kstage,
              Var<double,NUMCOMPS>& a_U,
              int a_stage
              )
{
  for (int n = 0; n < NUMCOMPS;n++)
    {
      if (a_stage == 0) 
        {
          a_UTaylor(n*RKORDER) = a_U(n);
          a_UTaylor(1 + n*RKORDER) = a_kstage(n);
          a_UTaylor(2 + n*RKORDER) = -(3.0/2)*a_kstage(n);
          a_UTaylor(3 + n*RKORDER) = (2.0/3)*a_kstage(n);
        }
      if (a_stage == 1)
        {
          a_UTaylor(2 + n*RKORDER) += a_kstage(n);
          a_UTaylor(3 + n*RKORDER) += (-2.0/3)*a_kstage(n);
        }
      if (a_stage == 2)
        {
          a_UTaylor(2 + n*RKORDER) += a_kstage(n);
          a_UTaylor(3 + n*RKORDER) += (-2.0/3)*a_kstage(n);
          
        }
      if (a_stage == 3)
        {
          a_UTaylor(2 + n*RKORDER) += (-1.0/2)*a_kstage(n);
          a_UTaylor(3 + n*RKORDER) += (2.0/3)*a_kstage(n);
        }
    }
}
PROTO_KERNEL_END(RKTaylorCoefsF, RKTaylorCoefs)
PROTO_KERNEL_START
void
RKTimeInterpF(Var<double,NUMCOMPS>& a_U,
              Var<double,RKORDER*NUMCOMPS>& a_UTaylor,
              double& a_chi
              )
{
  for (int n = 0; n < NUMCOMPS;n++)
    {
      a_U(n) = 0.;
      for (int term = 0; term < RKORDER; term++)
        {
          a_U(n) = a_U(n)*a_chi + a_UTaylor(RKORDER-term-1 + n*RKORDER);
        }
    }
}
PROTO_KERNEL_END(RKTimeInterpF,RKTimeInterp)
PROTO_KERNEL_START
void stageUpdateF(
                 Var<double,NUMCOMPS>& a_UStage,
                 Var<double,NUMCOMPS>& a_rhsTot,
                 const Var<double,NUMCOMPS>& a_U0,
                 const Var<double,NUMCOMPS>& a_kStage,
                 const double& a_dtStageNext,
                 const double& a_stageWeight)
{
  for (int comp = 0; comp < NUMCOMPS; comp++)
    {
      a_rhsTot(comp) += a_stageWeight*a_kStage(comp);
      a_UStage(comp) = a_U0(comp) + a_dtStageNext*a_kStage(comp); 
    }
}
PROTO_KERNEL_END(stageUpdateF,stageUpdate)
template
<class OpType, typename T,unsigned int C,MemType MEM>
AMRSubcycleExplicit<OpType,T,C,MEM>::AMRSubcycleExplicit
                                    (shared_ptr<AMRData<T,C,MEM> > a_dataPtr,
                                     const T& a_dx,
                                     int a_timeRefRatio,
                                     int a_level
                                     )
{
  define(a_dataPtr,a_dx,a_timeRefRatio,a_level);
}
template
<class OpType, typename T,unsigned int C,MemType MEM>
void AMRSubcycleExplicit<OpType,T,C,MEM>::define
                                   (shared_ptr<AMRData<T,C,MEM> > a_dataPtr,
                                    const T& a_dx,
                                    int a_timeRefRatio,
                                    int a_level
                                    )
{
  // Define data members for any level.
  //cout << "defining level " << a_level << endl;
  m_dataPtr = a_dataPtr;
  m_level = a_level;
  m_dbl = (*m_dataPtr)[m_level].layout();
  m_time = 0.;
  m_dx = a_dx;
  m_numSteps = 0;

  m_ghost = m_op.ghostSize();
  //cout << "ghostsize = " << endl;
  
  // If we have a finer level, make recursive call to constructor, set up storage
  // for the saving the Taylor expansion coefficients for the current level.
  
  if (m_level < (*m_dataPtr).numLevels() - 1)
    {
      T dxFine = m_dx/PR_AMR_REFRATIO;
      m_refRatioFine = PR_AMR_REFRATIO*Point::Ones();
      m_timeRefRatioFine = PR_AMR_REFRATIO;
      cout << "ref ratio fine" << m_timeRefRatioFine << endl ;
      m_finePtr = shared_ptr<AMRSubcycleExplicit<OpType,T,C,MEM> >
                            (new AMRSubcycleExplicit<OpType,T,C,MEM>
                             (m_dataPtr,
                              dxFine,
                              m_timeRefRatioFine,
                              m_level + 1)
                             );
      
      m_UTaylor.define(m_dbl,Point::Zeros());
      m_dblFine = (*m_dataPtr)[m_level+1].layout();
      m_register.define(m_dbl,m_dblFine,m_refRatioFine);
    }
  // If we have a coarser level:
  // (1) build interpolation stencil,
  // (2) set up storage for Taylor expansion coefficients from the next coarser level,
  // (3) For each patch, compute boxes over which the ghost cell interpolation will
  // be applied.
  if (m_level > 0)
    {
      m_dblCoarse = (*m_dataPtr)[m_level-1].layout();
      m_refRatioCoarse = PR_AMR_REFRATIO*Point::Ones();
      m_timeRefRatioCoarse = PR_AMR_REFRATIO;
      cout << "ref ratio coarse" << m_timeRefRatioCoarse << endl ;
      // Set up for fourth-order in space and time c/f bcs.

      // Storage for RK4 interpolation function
      Point ghostInterp = m_op.ghostSize()/m_refRatioCoarse + (4*Point::Ones());
      DisjointBoxLayout dblCoarsened = m_dbl.coarsen(m_refRatioCoarse);
      m_UTaylorCoarsened.define(
                                dblCoarsened,
                                ghostInterp);
      
      // Construct fourth-order spatial stencil.
      m_cfInterp = InterpStencil<T>::Build(3,Box(Point::Ones(-2),Point::Ones(2)),3,PR_AMR_REFRATIO);
    }
}
template
<class OpType, typename T,unsigned int C,MemType MEM>
void AMRSubcycleExplicit<OpType,T,C,MEM>::advance
(LevelFluxRegister<T,C,MEM>& a_coarseRegister,
                               const T& a_dt,
                               const T& a_chi,
                               bool a_willRegridCoarse)
{
  // cout << "entering advance, time = " << m_time << endl;
  // Setup for interpolation from coarse grid data using dense representation.
  PR_TIMERS("RK Advance");
  LevelBoxData<T,C,MEM> UStage(m_dbl,m_ghost);
  LevelBoxData<T,C,MEM> rhsTot(m_dbl,Point::Zeros());
  rhsTot.setToZero();
  (*m_dataPtr)[m_level].copyTo(UStage);
  LevelBoxData<T,C,MEM>& U0 = (*m_dataPtr)[m_level];
  if (m_level <  (*m_dataPtr).numLevels()-1) m_register.reset();
  for (int stage = 0; stage < RKSTAGES; stage++) // Stage loop.
    {     
      // Interpolate ghost cells. 
      if (m_level > 0)
        {
          T chistage = a_chi + m_stageTimes[stage]/m_timeRefRatioCoarse;
          T timeStage = m_time + a_dt*m_stageTimes[stage];
          for (auto dit=m_dbl.begin();*dit!=dit.end();++dit)
            {
              PR_TIMERS("Coarse-fine interpolation");
              BoxData<T,C,MEM> UStageBoxTemp(m_dbl[*dit]);
              auto & UStageBox = UStage[*dit];
              UStageBox.copyTo(UStageBoxTemp);             
              auto UTimeInterp = forall<T,C>
                (RKTimeInterp,m_UTaylorCoarsened[*dit],chistage);
              UStageBox |= m_cfInterp(UTimeInterp);
              UStageBoxTemp.copyTo(UStageBox,m_dbl[*dit]);
            }
        }
      UStage.exchange();
      T dtStageNext = a_dt*m_stageTimes[stage+1];
      T dtStage = a_dt*m_stageTimes[stage];
      
      for (auto dit = m_dbl.begin();*dit != dit.end();++dit) 
        {
          // interpolate coarse grid data in time, then in space, to fill ghost cells.
          auto & UStageBox = UStage[*dit];
          array<BoxData<T,C,MEM>, DIM> fluxes;
          
          for (int dir = 0; dir < DIM; dir++)
            { 
              fluxes[dir].define(m_dbl[*dit].extrude(dir));
            }

          BoxData<T,C,MEM> kStage(m_dbl[*dit]);     
          
          // Perform stage calculation and recompute UStage.
          {
            PR_TIMERS("RK Stage Update");
            m_op(UStageBox,kStage,fluxes,m_dx);
            
            forallInPlace(stageUpdate,
                          UStageBox,
                          rhsTot[*dit],
                          U0[*dit],
                          kStage,
                          dtStageNext,
                          m_stageWgt[stage]);
          }
          
          if (m_level < (*m_dataPtr).numLevels()-1)
            {
              // Store stage information.
              kStage *= a_dt;
              
              forallInPlace(RKTaylorCoefs,m_UTaylor[*dit],kStage,U0[*dit],stage);
              
              // Increment flux registers.
              
              T dtwgt = a_dt*m_stageWgt[stage];           
              for (int dir = 0; dir < DIM; dir++)
                {
                  m_register.incrementCoarse(fluxes[dir],*dit,dtwgt,dir);
                }
            }
          if (m_level> 0)
            {
              // Increment fine flux registers.
              
              T dtwgt = a_dt*m_stageWgt[stage];
              for (int dir = 0; dir < DIM; dir++)
                {
                  a_coarseRegister.incrementFine(fluxes[dir],*dit,dtwgt,dir);
                }
            }
        } // end loop for a single stage.
     
    } // end loop over RK4 stages.

  for (auto dit = m_dbl.begin();*dit != dit.end();++dit)
    {
      PR_TIMERS("RK Update");
      rhsTot[*dit] *= a_dt;
      U0[*dit] += rhsTot[*dit];
    }

  // recursive call to advance finer levels.
  if (m_level < (*m_dataPtr).numLevels() - 1)
    {
  // Copy Taylor coefs to the data holder on the next finer level.
      m_UTaylor.copyTo(m_finePtr->m_UTaylorCoarsened);

  // Advance finer levels.
      m_numSteps++;
      bool willRegrid = (m_numSteps == m_regridInterval - 1);
      // cout << "loop for next finer " << endl;
      for (int nstep = 0; nstep < m_timeRefRatioFine;nstep++)
        {
          // cout << "fine nstep = " << nstep << endl;
          double chi = nstep*1.0/m_timeRefRatioFine;
          double dtFine = a_dt/m_timeRefRatioFine;
          m_finePtr->advance(m_register,
                             dtFine,
                             chi,
                             willRegrid);
        }
  
  // Average down and reflux.
      double dxI = 1.0/m_dx;
      LevelBoxData<T,C,MEM> UCoarsened(m_dblFine.coarsen(m_refRatioFine),Point::Zeros());
      Stencil<double> avgdown =Stencil<double>::AvgDown(m_refRatioFine);
      
      for (auto dit = UCoarsened.begin();*dit != dit.end(); ++ dit)
        {
          UCoarsened[*dit] |= avgdown((*m_dataPtr)[m_level+1][*dit]);
        }
      UCoarsened.copyTo((*m_dataPtr)[m_level]);
      // Minus sign in the coefficient since the RHS is -div(F).
      m_register.reflux((*m_dataPtr)[m_level],-dxI);
    }
  
  m_time += a_dt;
  
  // Regrid. Regridding will be managed by the finest level that doesn't change.
  if ((m_level < (*m_dataPtr).numLevels() - 1)&&(m_numSteps == m_regridInterval-1))
    {
      m_numSteps = 0;
      if (!a_willRegridCoarse)
        {
          // call regridding - TBD
          cout << "Warning - you've hit regridding which is not yet implemented" << endl;
        }
    }
}
