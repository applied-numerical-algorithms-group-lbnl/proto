
PROTO_KERNEL_START
void
f_RKTaylorCoefsF(Var<double,RKORDER*(NUMCOMPS)>& a_UTaylor,
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
PROTO_KERNEL_END(f_RKTaylorCoefsF, f_RKTaylorCoefs)

PROTO_KERNEL_START
void
f_RKTimeInterpF(Var<double,NUMCOMPS>& a_U,
                Var<double,RKORDER*(NUMCOMPS)>& a_UTaylor,
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
PROTO_KERNEL_END(f_RKTimeInterpF,f_RKTimeInterp)

PROTO_KERNEL_START
void f_stageUpdateF(
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
PROTO_KERNEL_END(f_stageUpdateF,f_stageUpdate)

template <class OpType, typename T,unsigned int C,MemType MEM>
AMRRK4<OpType,T,C,MEM>::AMRRK4(
        AMRData<T, C, MEM>& a_data,
        const Array<T, DIM>& a_dx,
        int a_regridInterval,
        int a_bufferSize)
{
  define(a_data, a_dx, 0,a_regridInterval,a_bufferSize);
}

template<typename T>
void trimVector(std::vector<T>& a_vector, unsigned int a_size)
{
    while (a_vector.size() > a_size)
    {
        a_vector.pop_back();
    }
}

template <class OpType, typename T,unsigned int C,MemType MEM>
void AMRRK4<OpType,T,C,MEM>::define(
                                    AMRData<T, C, MEM>& a_data,
                                    const Array<T, DIM>& a_dx,
                                    int a_regridInterval,
                                    int a_bufferSize)
{
  define(a_data, a_dx, 0 ,a_regridInterval,a_bufferSize);
}

template <class OpType, typename T,unsigned int C,MemType MEM>
void AMRRK4<OpType,T,C,MEM>::define(
                                    AMRData<T, C, MEM>& a_data,
                                    const Array<T, DIM>& a_dx,
                                    int                 a_level,
                                    int                 a_regridInterval,
                                    int                 a_bufferSize)
{
     PR_TIMERS("AMRRK4::define");
    // dont do anything fancy, just kill everything and redefine
    if (a_level == 0)
    {
        int numLevels = a_data.grid().numLevels();
        
        m_UStage.clear();
        m_RHSTotal.clear();
        m_UTaylor.clear();
        m_UTaylorC.clear();
        m_fluxRegisters.clear();
        m_ops.clear();
        
        m_UStage.resize(numLevels);
        m_RHSTotal.resize(numLevels);
        m_UTaylor.resize(numLevels - 1);
        m_UTaylorC.resize(numLevels - 1);
        m_fluxRegisters.resize(numLevels - 1);       

        m_ops.resize(numLevels);
        if (m_numSteps.size() != numLevels)
          {
            m_numSteps= std::vector<int>(numLevels,0);
          }
        m_timeRefRatios.resize(numLevels - 1);
        m_spaceRefRatios.resize(numLevels - 1);
        m_regridInterval.resize(numLevels - 1);
        m_bufferSize.resize(numLevels - 1);
        m_cfInterp.resize(numLevels - 1);
        
        m_data = &a_data;
        Array<T, DIM> dx = a_dx;
        for (int lvl = a_level; lvl < numLevels; lvl++)
        {
            auto& layout = a_data.grid()[lvl];

            m_ops[lvl].define(dx);
            m_UStage[lvl] = std::make_shared<LevelBoxData<T, C, MEM>>(layout, OpType::ghost());
            m_UStage[lvl]->setToZero();
            m_RHSTotal[lvl] = std::make_shared<LevelBoxData<T, C, MEM>>(layout, Point::Zeros());
            m_RHSTotal[lvl]->setToZero();
            if (lvl < numLevels-1)
            {
                auto& fineLayout = a_data.grid()[lvl+1];
                Point spaceRefRatio = a_data.grid().refRatio(lvl);
                auto cfLayout = fineLayout.coarsen(spaceRefRatio);

                m_regridInterval[lvl] = a_regridInterval;
                m_bufferSize[lvl] = a_bufferSize;
                m_spaceRefRatios[lvl] = spaceRefRatio;
                m_timeRefRatios[lvl] = spaceRefRatio.max();
                m_cfInterp[lvl] = InterpStencil<T>::FiniteVolume(spaceRefRatio, OpType::order());
                m_UTaylor[lvl] =
                    std::make_shared<LevelBoxData<T, C*RKORDER, MEM>>(layout, Point::Zeros());
                m_UTaylor[lvl]->setToZero();
                m_fluxRegisters[lvl] =
                    std::make_shared<LevelFluxRegister<T, C, MEM>>(layout, fineLayout, spaceRefRatio, dx);

                m_fluxRegisters[lvl]->reset();
                // Point::Ones() are here to make the division a ceil instead of a floor
                Point interpGhost = OpType::ghost() / spaceRefRatio + Point::Ones();
                interpGhost += m_cfInterp[lvl].ghost();

                m_UTaylorC[lvl] =
                    std::make_shared<LevelBoxData<T, C*RKORDER, MEM>>(cfLayout, interpGhost);
                m_UTaylorC[lvl]->setToZero();
                for (int dir = 0; dir < DIM; dir++)
                {
                    dx[dir] /= spaceRefRatio[dir];
                }
            }
        }
    }
    // leave coarse level stuff alone and only redefine things on the finer levels
    else 
    {
        int numLevels = a_data.grid().numLevels();
        trimVector(m_UStage,        a_level);
        trimVector(m_RHSTotal,      a_level);
        trimVector(m_UTaylor,       a_level);
        trimVector(m_UTaylorC,      a_level);
        trimVector(m_fluxRegisters, a_level);

        m_ops.resize(numLevels);
        //m_numSteps.resize(numLevels);
        m_regridInterval.resize(numLevels-1);
        m_bufferSize.resize(numLevels-1);
        m_timeRefRatios.resize(numLevels - 1);
        m_spaceRefRatios.resize(numLevels - 1);
        m_cfInterp.resize(numLevels - 1);
        
        m_data = &a_data;

        Array<T, DIM> dx = a_dx;
        for (int lvl = a_level; lvl < numLevels; lvl++)
        {
            auto& layout = a_data.grid()[lvl];
            m_ops[lvl].define(dx);
            m_UStage.push_back(std::make_shared<LevelBoxData<T, C, MEM>>(layout, OpType::ghost()));
            m_RHSTotal.push_back(std::make_shared<LevelBoxData<T, C, MEM>>(layout, Point::Zeros()));
            //m_numSteps[lvl] = 0;
            if (lvl < numLevels-1)
            {
                auto& fineLayout = a_data.grid()[lvl+1];
                Point spaceRefRatio = a_data.grid().refRatio(lvl);
                auto cfLayout = fineLayout.coarsen(spaceRefRatio);

                m_spaceRefRatios[lvl] = spaceRefRatio;
                m_timeRefRatios[lvl] = spaceRefRatio.max();
                m_regridInterval[lvl] = a_regridInterval;
                m_bufferSize[lvl] = a_bufferSize;
                m_cfInterp[lvl] = InterpStencil<T>::FiniteVolume(spaceRefRatio, OpType::order());
                m_UTaylor.push_back(
                        std::make_shared<LevelBoxData<T, C*RKORDER, MEM>>(layout, Point::Zeros()));
                m_fluxRegisters.push_back(
                        std::make_shared<LevelFluxRegister<T, C, MEM>>(layout, fineLayout, spaceRefRatio, dx));

                // Point::Ones() are here to make the division a ceil instead of a floor
                Point interpGhost = OpType::ghost() / spaceRefRatio + Point::Ones();
                interpGhost += m_cfInterp[lvl].ghost();

                m_UTaylorC.push_back(
                        std::make_shared<LevelBoxData<T, C*RKORDER, MEM>>(cfLayout, interpGhost));
                for (int dir = 0; dir < DIM; dir++)
                {
                    dx[dir] /= spaceRefRatio[dir];
                }
            }
        }
        PROTO_ASSERT(m_UStage.size() == numLevels, 
                "AMRRK4::define | Error: Wrong number of elements.")
            PROTO_ASSERT(m_RHSTotal.size() == numLevels, 
                    "AMRRK4::define | Error: Wrong number of elements.")
            PROTO_ASSERT(m_UTaylor.size() == numLevels - 1, 
                    "AMRRK4::define | Error: Wrong number of elements.")
            PROTO_ASSERT(m_UTaylorC.size() == numLevels - 1, 
                    "AMRRK4::define | Error: Wrong number of elements.")
            PROTO_ASSERT(m_fluxRegisters.size() == numLevels - 1, 
                    "AMRRK4::define | Error: Wrong number of elements.")
    }

}

template <class OpType, typename T,unsigned int C,MemType MEM>
void AMRRK4<OpType,T,C,MEM>::advance(T a_dt)
{
    advance(0, a_dt, 0);
}

template <class OpType, typename T,unsigned int C,MemType MEM>
void AMRRK4<OpType,T,C,MEM>::advance(int a_level, T a_dt, T a_chi)
{
#ifdef PR_HDF5
  HDF5Handler h5;
#endif
    PR_TIMERS("AMRRK4::advance");
    auto dx = m_ops[a_level].dx();
    auto dx0 = m_ops[0].dx();

    auto& layout    = m_data->grid()[a_level];
    auto& U0 = (*m_data)[a_level];
    auto& UStage = *m_UStage[a_level];
    auto& RHSTotal = *m_RHSTotal[a_level];
    
    

    U0.copyTo(UStage);
    RHSTotal.setToZero();
    
    for (int stage = 0; stage < RKSTAGES; stage++)
    {
        if (a_level > 0)
        {       
          interpBoundaries(a_level, stage, a_chi);
        }
        
        stageUpdate(a_level, stage, a_dt);
       
    } // end stage loop
    
    // Update U0
    for (auto iter : layout)
    {
        RHSTotal[iter] *= a_dt;
        U0[iter] += RHSTotal[iter];
    }

    // Recursive Advance
    if (a_level < m_data->numLevels() - 1)
    {
      int timeRefRatio = m_timeRefRatios[a_level];
      // Copy Taylor coeffs in time to coarse grid data holder aligned with finer level.
      m_UTaylor[a_level]->copyTo(*m_UTaylorC[a_level]);
      for (int nstep = 0; nstep < timeRefRatio; nstep++)
        {
          T chi = nstep*1.0/timeRefRatio;
          T dtFine = a_dt/timeRefRatio;
          advance(a_level+1, dtFine, chi);
        }
      Proto::averageDown((*m_data)[a_level], (*m_data)[a_level+1], m_spaceRefRatios[a_level]);
      // 1/dx factor is done at the increment stage to account for anisotropy
      {
        PR_TIMERS("AMRRK4::advance reflux");
        m_fluxRegisters[a_level]->reflux(U0, 1.0);
        m_fluxRegisters[a_level]->reset();
      }
    }
    m_numSteps[a_level]++;
    // if (a_level == 1) cout << m_numSteps[a_level] << endl;
    // Do Regridding
    if (a_level < m_data->numLevels() - 1)
      {
        bool regridThisLevel = false;
        if (Proto::procID()==0)
          {
            //std::cout <<"regrid level: " << a_level << " , " << m_numSteps[a_level] <<
            //  " , " << m_regridInterval[a_level] << std::endl;
          }
        if (a_level  == 0)
          {
            regridThisLevel = (m_numSteps[0]%m_regridInterval[0] == 0);
          }
        else
          {                              
            regridThisLevel = (m_numSteps[a_level]%m_regridInterval[a_level] == 0) 
              && ((m_numSteps[a_level])%
                  ((m_timeRefRatios[a_level-1])*m_regridInterval[a_level-1]) != 0);
          }  
        if (regridThisLevel)
          {
            if (Proto::procID()==0)
              {
                //std::cout <<" regridding: " << a_level << ", step = " << m_numSteps[a_level] <<  std::endl;
              }
            //h5.writeAMRData(dx0,*m_data,
            //                      "U_lev_pre"+to_string(a_level)+"_NC_"+to_string(m_data->grid()[0].domain().sizes()[0])+ "_step_"+to_string(m_numSteps[a_level]));
            
            regrid(a_level);
            // h5.writeAMRData(dx0,*m_data,
            //              "U_lev_post"+to_string(a_level)+"_NC_"+to_string(m_data->grid()[0].domain().sizes()[0])+ "_step_"+to_string(m_numSteps[a_level]));
       
          }
      }
    else
      {
        //h5.writeAMRData(dx0,*m_data,
        //                    "U_lev_"+to_string(a_level)+"_NC_"+to_string(m_data->grid()[0].domain().sizes()[0])+ "_step_"+to_string(m_numSteps[a_level]));
      }
}

template <class OpType, typename T,unsigned int C,MemType MEM>
void AMRRK4<OpType,T,C,MEM>::interpBoundaries(int a_level, int a_stage, T a_chi)
{
    PROTO_ASSERT(a_level > 0,
        "AMRRK4::interpBoundaries | Error: Interpolation is only done for level > 0.");
   
#ifdef PR_HDF5
    //HDF5Handler h5;
#endif
     PR_TIMERS("AMRRK4::interpBoundaries");
    auto dx = m_ops[a_level].dx();

    auto& layout    = m_data->grid()[a_level];
    auto& UStage    = *m_UStage[a_level];
    auto& UTaylorC  = *m_UTaylorC[a_level-1];
    auto& cfInterp  = m_cfInterp[a_level-1];

    T chi = a_chi + m_stageTimes[a_stage]/m_timeRefRatios[a_level-1];
    Stencil<T> I = 1.0*Shift::Zeros();

    for (auto iter : layout)
    {
        if (layout.onLevelBoundary(layout.point(iter)))
        {
            // get data
            auto& UStage_i   = UStage[iter];
            auto& UTaylorC_i = UTaylorC[iter];
            BoxData<T, C, MEM> UStageTemp(layout[iter]);
            
            // save copy of interior
            UStageTemp |= I(UStage_i);

            // interpolate in time
            auto UTimeInterp = forall<T, C>(f_RKTimeInterp, UTaylorC_i, chi);
        
            // interpolate boundary conditions
            {
              PR_TIMERS("InterpBoundaries single patch");
              UStage_i |= cfInterp(UTimeInterp);
            }
            //cfInterp.apply(UStage_i, UTimeInterp);
        
            // copy back the interior
            UStage_i |= I(UStageTemp);
        }
    }
    UStage.exchange();

}

template <class OpType, typename T,unsigned int C,MemType MEM>
void AMRRK4<OpType,T,C,MEM>::stageUpdate(int a_level, int a_stage, T a_dt)
{
#ifdef PR_HDF5
  //HDF5Handler h5;
#endif
     PR_TIMERS("AMRRK4::stageUpdate");
    
    auto dx = m_ops[a_level].dx();

    auto& layout    = m_data->grid()[a_level];
    auto& UStage    = *m_UStage[a_level];
    auto& U0        = (*m_data)[a_level];
    auto& RHSTotal  = *m_RHSTotal[a_level];
    auto& op        = m_ops[a_level];
    
    T dtStageNext = a_dt*m_stageTimes[a_stage+1];
    T stageWeight = m_stageWeights[a_stage];
    T dtWeight    = a_dt*m_stageWeights[a_stage];   
    op.bcStage(UStage,U0,a_stage);
    
    for (auto iter : layout)
    {
        // get data
        auto& UStage_i   = UStage[iter];
        auto& U0_i       = U0[iter];
        auto& RHSTotal_i = RHSTotal[iter];

        Array<BoxData<T, C, MEM>, DIM> fluxes;
        for (int dir = 0; dir < DIM; dir++)
        {
            fluxes[dir].define(layout[iter].grow(dir, Side::Hi, 1));
            fluxes[dir].setToZero();
        }
        BoxData<T,C,MEM> kStage(layout[iter], 0.0);
        
        // apply spatial operator to compute k
        {
          PR_TIMERS("AMRRK4::stageUpdate apply operator");
          op(kStage, fluxes, UStage_i);
        }
        for (int dir = 0; dir < DIM ;dir++)
        {
          //BoxData<T> flux_d = slice(Fluxes[iter], dir);
          //fluxes[dir].copyTo(flux_d);
        }

        // do stage update
        {
           PR_TIMERS("AMRRK4::stageUpdate do stage update");
           forallInPlace(f_stageUpdate, UStage_i, RHSTotal_i, U0_i, kStage, dtStageNext, stageWeight);
        }
        if (a_level < m_data->numLevels() - 1)
        {
            auto& UTaylor_i  = (*m_UTaylor[a_level])[iter];
            // store stage information
            kStage *= a_dt;
            forallInPlace(f_RKTaylorCoefs, UTaylor_i, kStage, U0_i, a_stage);

            // increment flux register
            for (int dir = 0; dir < DIM; dir++)
            {
                T cdx = op.dx()[dir];
                m_fluxRegisters[a_level]->incrementCoarse(fluxes[dir], iter, dtWeight, dir);
            }
        }
        if (a_level > 0)
        {
            for (int dir = 0; dir < DIM; dir++)
            {
                T cdx = m_ops[a_level-1].dx()[dir];
                m_fluxRegisters[a_level-1]->incrementFine(fluxes[dir], iter, dtWeight, dir);
            }
        }
    } // end patch loop
}

template <class OpType, typename T,unsigned int C,MemType MEM>
void AMRRK4<OpType,T,C,MEM>::regrid(int a_baseLevel)
{
#ifdef PR_HDF5
  HDF5Handler h5;
#endif
    PR_TIMERS("AMRRK4::regrid");
    AMRGrid newGrids = m_data->grid();

    // for all coarse levels (finest to coarsest)
    for (int level = m_data->numLevels()-2; level >= a_baseLevel; level--)
    {
        auto& layout = newGrids[level];

        // Compute tags on this level
        LevelTagData tags(layout, Point::Ones(m_bufferSize[level]));
        tags.setToZero();
        for (auto iter : layout)
        {
            auto& tags_i = tags[iter];
            auto& U0_i   = (*m_data)[level][iter];
            m_ops[level].generateTags(tags_i, U0_i);
        }
        tags.setBoundary(0);
        newGrids.addFinerTags(tags, level);
        AMRGrid::buffer(tags, Point::Ones(m_bufferSize[level]));
        newGrids.regrid(tags, level);
    }
    for (int level = max(2, a_baseLevel); level < newGrids.numLevels(); level++)
    {
        newGrids.enforceNesting2(level);
    }
    auto dx0 = m_ops[0].dx();
     int a_level = a_baseLevel;
     //h5.writeAMRData(dx0,*m_data,
     //                        "U_lev_definePre"+to_string(a_level)+"_NC_"+to_string(m_data->grid()[0].domain().sizes()[0])+ "_step_"+to_string(m_numSteps[a_level]));
    
    m_data->regrid(newGrids, a_baseLevel, OpType::order());
    
    //h5.writeAMRData(dx0,*m_data,
    //                         "U_lev_definePost"+to_string(a_level)+"_NC_"+to_string(m_data->grid()[0].domain().sizes()[0])+ "_step_"+to_string(m_numSteps[a_level]));
    
    define(*m_data, m_ops[a_baseLevel].dx(), a_baseLevel,
           m_regridInterval[a_baseLevel],m_bufferSize[a_baseLevel]);            
}
