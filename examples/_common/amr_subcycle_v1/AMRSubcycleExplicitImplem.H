inline void writeAMRGridPlot2
(double a_dx,AMRGrid a_grids,string a_label)
  {
    AMRData<short,1> data(a_grids,Point::Zeros());
    for (int ii = 0; ii < a_grids.numLevels();ii++)
      {
        auto & dataLev = data[ii];
        for (auto dit=dataLev.begin();*dit != dit.end();++dit)
          {
            dataLev[*dit].setVal(ii);
          }
      }
    HDF5Handler h5;
    h5.writeAMRData(a_dx,data,a_label);
  }
inline void
writeFluxRegisterContents(LevelFluxRegister<double,NUMCOMPS>& a_reg,
                           string a_str,double a_dx)
{
  DisjointBoxLayout dblc = a_reg.crseLayout();
  LevelBoxData<double,NUMCOMPS> test(dblc,Point::Zeros());
  test.setToZero();
  a_reg.reflux(test,1.0);
  HDF5Handler h5;
  h5.writeLevel(a_dx,test,a_str);
}
void
writeFluxRegisterForm(LevelFluxRegister<double,NUMCOMPS>& a_reg,
                           string a_str,double a_dx)
{
  a_reg.reset(1.);
  writeFluxRegisterContents(a_reg,a_str,a_dx);
  a_reg.reset();
}                 
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
                                     shared_ptr<vector<int > > a_globalStep,
                                     const T& a_dx,
                                     int a_timeRefRatio,
                                     int a_level
                                     )
{
  define(a_dataPtr, a_globalStep, a_dx, a_timeRefRatio, a_level);
}
template
<class OpType, typename T,unsigned int C,MemType MEM>
void AMRSubcycleExplicit<OpType,T,C,MEM>::define
                                   (shared_ptr<AMRData<T,C,MEM> > a_dataPtr,
                                    shared_ptr<vector<int > > a_globalStep,
                                    const T& a_dx,
                                    int a_timeRefRatio,
                                    int a_level
                                    )
{
  // Define data members for any level.
  //cout << "defining level " << a_level << endl;
  m_dataPtr = a_dataPtr;
  m_level = a_level;
  m_time = 0.;
  m_dx = a_dx;
  m_regridInterval = 1;
  m_globalStep = a_globalStep;
  m_ghost = m_op.ghostSize();
  
  // If we have a finer level, make recursive call to constructor, set up storage
  // for the saving the Taylor expansion coefficients for the current level.
  
  if (m_level < (*m_dataPtr).numLevels() - 1)
    {
      auto dbl = m_dataPtr->grid()[m_level];
      // FIXME: Assumes isotropic refinement
      int spaceRefRatio = m_dataPtr->grid().refRatio(m_level)[0];
      T dxFine = m_dx/spaceRefRatio;
      m_refRatioFine = spaceRefRatio*Point::Ones();
      m_timeRefRatioFine = spaceRefRatio;
      m_finePtr = shared_ptr<AMRSubcycleExplicit<OpType,T,C,MEM> >
                            (new AMRSubcycleExplicit<OpType,T,C,MEM>
                             (m_dataPtr,
                              m_globalStep,
                              dxFine,
                              m_timeRefRatioFine,
                              m_level + 1)
                             );
      
      m_UTaylor.define(dbl,Point::Zeros());
      auto dblFine = (*m_dataPtr)[m_level+1].layout();
      PR_assert(dblFine == m_dataPtr->grid()[m_level+1]);
      m_register = shared_ptr<LevelFluxRegister<T,C,MEM> >
        (new LevelFluxRegister<T,C,MEM>(dbl,dblFine,m_refRatioFine));
    }
  // If we have a coarser level:
  // (1) build interpolation stencil,
  // (2) set up storage for Taylor expansion coefficients from the next coarser level,
  // (3) For each patch, compute boxes over which the ghost cell interpolation will
  // be applied.
  if (m_level > 0)
    {
      // FIXME: Assumes isotropic refinement
      int spaceRefRatio = m_dataPtr->grid().refRatio(m_level-1)[0];
      auto dbl = m_dataPtr->grid()[m_level];
      auto dblCoarse = m_dataPtr->grid()[m_level-1];
      m_refRatioCoarse = spaceRefRatio*Point::Ones();
      m_timeRefRatioCoarse = spaceRefRatio;
      // cout << "ref ratio coarse" << m_timeRefRatioCoarse << endl ;
      // Set up for fourth-order in space and time c/f bcs.

      // Storage for RK4 interpolation function
      Point ghostInterp = m_op.ghostSize()/m_refRatioCoarse + (4*Point::Ones());
      DisjointBoxLayout dblCoarsened = dbl.coarsen(m_refRatioCoarse);
      m_UTaylorCoarsened.define(
                                dblCoarsened,
                                ghostInterp);
      //cout << "in define: level = " << m_level << endl;
      //m_UTaylorCoarsened.layout().print();
      // Construct fourth-order spatial stencil.
      m_cfInterp = InterpStencil<T>::Build(4, spaceRefRatio);
      m_cfInterpTP.define(5,spaceRefRatio);
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
  
  auto dbl = m_dataPtr->grid()[m_level];
  PR_assert(dbl == (*m_dataPtr)[m_level].layout());
  int ind;
  DisjointBoxLayout dblFiner,dblCoarser;
  if (m_level < m_dataPtr->numLevels()-1)
    {
      dblFiner = m_dataPtr->grid()[m_level+1];
      PR_assert(dblFiner == (*m_dataPtr)[m_level+1].layout());
    }
  if (m_level > 0)
    {
      dblCoarser = m_dataPtr->grid()[m_level-1];
      PR_assert(dblCoarser == (*m_dataPtr)[m_level-1].layout());
    }
  LevelBoxData<T,C,MEM> UStage(dbl,m_ghost);
  (*m_dataPtr)[m_level].copyTo(UStage);
  LevelBoxData<T,C,MEM> rhsTot(dbl,Point::Zeros());
  LevelBoxData<T,C,MEM>& U0 = (*m_dataPtr)[m_level];
  PR_assert((U0.layout() == dbl));
  if (m_level <  (*m_dataPtr).numLevels()-1) m_register->reset();
  for (int stage = 0; stage < RKSTAGES; stage++) // Stage loop.
    {     
      // Interpolate ghost cells.
      if (m_level > 0)
        {
          T chistage = a_chi + m_stageTimes[stage]/m_timeRefRatioCoarse;
          T timeStage = m_time + a_dt*m_stageTimes[stage];
          int counter = 0;
          for (auto dit=dbl.begin();*dit!=dit.end();++dit)
            {
              //Shutting off interpolation if this patch issurrounded by other fine patches.
              auto ptcen = dbl.point(*dit);
              bool skipInterp = true;
              Box bx(ptcen - Point::Ones(),ptcen + Point::Ones());
              for (auto bxit = bx.begin();bxit.ok();++bxit)
                {
                  if (dbl.find(*bxit) == dbl.end())
                    {
                      skipInterp = false;
                      break;
                    }
                }
              if (!skipInterp)
                {
                  PR_TIMERS("Coarse-fine interpolation");
                  
                  BoxData<T,C,MEM> UStageBoxTemp(dbl[*dit]);
                  Stencil<T> idOp = 1.0*Shift::Zeros();
                  auto & UStageBox = UStage[*dit];
                  
                  UStageBox.copyTo(UStageBoxTemp);
                  auto UTimeInterp = forall<T,C>
                    (RKTimeInterp,m_UTaylorCoarsened[*dit],chistage);
                  { 
                    PR_TIMERS("TP Coarse-fine interpolation in space");
                    m_cfInterpTP.apply(UStageBox,UTimeInterp);
                  }
                  
                  {
                    PR_TIMERS("final copy");
                    UStageBox |= idOp(UStageBoxTemp);                  
                  }
                }
            }
        }

      UStage.exchange();
      T dtStageNext = a_dt*m_stageTimes[stage+1];
      T dtStage = a_dt*m_stageTimes[stage];
      for (auto dit = dbl.begin();*dit != dit.end();++dit) 
        {
              
          auto & UStageBox = UStage[*dit];
          array<BoxData<T,C,MEM>, DIM> fluxes;
              
          for (int dir = 0; dir < DIM; dir++)
            { 
              fluxes[dir].define(dbl[*dit].extrude(dir));
              fluxes[dir].setVal(0.);
            }
              
          BoxData<T,C,MEM> kStage(dbl[*dit]);     
          kStage.setVal(0.);
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
              // Increment fine flux registers.              
              for (int dir = 0; dir < DIM; dir++)
                {
                  m_register->incrementCoarse(fluxes[dir],*dit,dtwgt,dir);
                }
            }
          if (m_level> 0)
            {
              // Increment coarse flux registers.
              T dtwgt = a_dt*m_stageWgt[stage];
              for (int dir = 0; dir < DIM; dir++)
                {  
                  a_coarseRegister.incrementFine(fluxes[dir],*dit,dtwgt,dir);
                }
            }
        } // end loop for a single stage.    
    } // end loop over RK4 stages.
  
  for (auto dit = dbl.begin();*dit != dit.end();++dit)
    {
      PR_TIMERS("RK Update");
      rhsTot[*dit] *= a_dt;
      U0[*dit] += rhsTot[*dit];      
    }
  // recursive call to advance finer levels
  bool willRegrid = ((m_numSteps+1)%m_regridInterval == 0);
  if (m_level + 1 < (*m_dataPtr).numLevels())
    {
      // Copy Taylor expansions coefficients to holder aligned
      // with next finer grids.
      
      m_UTaylor.copyTo(m_finePtr->m_UTaylorCoarsened);
        
      // Advance finer levels.
       
      for (int nstep = 0; nstep < m_timeRefRatioFine;nstep++)
        { 
          double chi = nstep*1.0/m_timeRefRatioFine;
          double dtFine = a_dt/m_timeRefRatioFine;
          bool willRegridLocal;
          if (nstep == m_timeRefRatioFine-1)
            {
              willRegridLocal = willRegrid;
            }
          else
            {
              willRegridLocal = false;
            }
          m_finePtr->advance(*m_register,
                             dtFine,
                             chi,
                             willRegridLocal);
        }
      // Average down and reflux.
      double dxI = 1.0/m_dx;
      LevelBoxData<T,C,MEM> UCoarsened(dblFiner.coarsen(m_refRatioFine),Point::Zeros());
      Stencil<double> avgdown =Stencil<double>::AvgDown(m_refRatioFine);
       
      for (auto dit = UCoarsened.begin();*dit != dit.end(); ++ dit)
        {
          UCoarsened[*dit] |= avgdown((*m_dataPtr)[m_level+1][*dit]);
        }
       
      UCoarsened.copyTo((*m_dataPtr)[m_level]);
       
      // Minus sign in the coefficient since the RHS is -div(F).
       
      m_register->reflux((*m_dataPtr)[m_level],-dxI);
    }
  m_time += a_dt;
  m_numSteps++;
  (*m_globalStep)[m_level]++;
  
  // Regrid. Regridding will be managed by the finest level that doesn't change.
  
  if ((m_level+1< (*m_dataPtr).numLevels())&&willRegrid)
    {
      if (!a_willRegridCoarse)
        {
         PR_TIMERS("Regrid"); 
          AMRGrid newgrids = (*m_dataPtr).grid();
          /* 
             starting at the finest level - 1
             (1) generate tags calling m_op function
             (2) union with tags at next finer level, if they exist.
             (3) generate grids at this level (w/o proper nesting) 
          */
          int buffersize=2;
          for (int lev = (*m_dataPtr).numLevels() - 2; lev >= m_level; lev--)
            {
              LevelTagData levtags(newgrids[lev],Point::Ones(buffersize));
              m_op.generateTags(levtags,(*m_dataPtr)[lev],m_dx,m_time);

              /*if ((lev == 1)&&(m_level == 0))
                {
                  HDF5Handler h5;
                  h5.writeLevel(m_dx/2,levtags,"Tags0PreBuffer%i",(*m_globalStep)[m_level]);
                  cout << "m_level, time step = " << m_level << " , " << (*m_globalStep)[m_level]<< endl;
                  }     */

              if ((lev + 2 <  (*m_dataPtr).numLevels())&&((*m_globalStep)[m_level] !=2 ))
                newgrids.addFinerGrids(levtags,lev);
              
              AMRGrid::buffer(levtags,Point::Ones(buffersize));

              /*if ((lev == 0)&&(m_level == 0))
                {
                  HDF5Handler h5;
                  h5.writeLevel(m_dx,levtags,"Tags00PostBuffer%i",(*m_globalStep)[m_level]);
                }     
              if ((lev == 1)&&(m_level == 1))
                {
                  HDF5Handler h5;
                  h5.writeLevel(m_dx,levtags,"Tags11PostBuffer%i",(*m_globalStep)[m_level]);
                  }*/
              newgrids.regrid(levtags,lev);
                       
            }
          if (m_level == 0)
            {
              //cout << "Level 0: numFinestGridsEnd = " << newgrids[2].size() << endl;
              //string str = "RegridLevel0PreNesting_"+to_string((*m_globalStep)[m_level]);
              //writeAMRGridPlot(m_dx,newgrids,str);
            }  
          for (int lev = max(2,m_level); lev < newgrids.numLevels(); lev++)
            {
              newgrids.enforceNesting2(lev);        
            }
          if (newgrids[2].size() == 0)
            {
              cout << "finest grids have disappeared !" << endl;
              cout << "m_level = " << m_level <<" , " << "step = "<< (*m_globalStep)[m_level] << endl;
              abort();
            }
          if (m_level == 0)
            {
              //cout << "Level 0: numFinestGridsEnd = " << newgrids[2].size() << endl;
              //string str = "RegridLevel0_"+to_string((*m_globalStep)[m_level]);
              //writeAMRGridPlot(m_dx,newgrids,str);
            }
          if (m_level == 1)
            {
              //cout << "Level 1: numFinestGrids = " << newgrids[2].size() << endl;
              //string str = "RegridLevel1_"+to_string((*m_globalStep)[m_level]);
              //writeAMRGridPlot(m_dx,newgrids,str);
            }
          // Regrid solution data for levels m_level+1, ... max level.
          
          m_dataPtr->regrid(newgrids,m_level,5);
          double dxLevel = m_dx;                       
          //Initialize data at all levels that have been regridded.
          this->define(m_dataPtr,m_globalStep,m_dx,m_timeRefRatioFine,m_level);
        }
    }
}
