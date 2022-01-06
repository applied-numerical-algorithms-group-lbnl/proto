
template <template<typename, MemType> class OPType, typename T, MemType MEM>
AMROp<OPType, T, MEM>::AMROp()
{
    m_defined = false;
}

template <template<typename, MemType> class OPType, typename T, MemType MEM>
AMROp<OPType, T, MEM>::AMROp(AMRGrid& a_grid, double a_cdx)
{
    define(a_grid, a_cdx);
}

template <template<typename, MemType> class OPType, typename T, MemType MEM>
void AMROp<OPType, T, MEM>::define(AMRGrid& a_grid, double a_cdx)
{
    int numLevels = a_grid.numLevels();
    m_dx.clear();
    m_levelOps.clear();
    m_fluxRegisters.clear();

    m_dx.resize(numLevels);
    m_levelOps.resize(numLevels);
    m_fluxRegisters.resize(numLevels - 1); //coarse levels only

    m_grid = a_grid; //intentional deep copy
    double dx = a_cdx;
    for (int lvl = 0; lvl < numLevels; lvl++)
    {
        m_dx[lvl] = dx;
        m_levelOps[lvl].define(dx);
        if (lvl < numLevels - 1)
        {
            m_fluxRegisters[lvl] =
                std::make_shared<LevelFluxRegister<T>>(
                    a_grid[lvl], a_grid[lvl+1],
                    Point::Ones(PR_AMR_REFRATIO)
                );
        }
        dx /= PR_AMR_REFRATIO;
    }
    m_defined = true;
}

template <template<typename, MemType> class OPType, typename T, MemType MEM>
void AMROp<OPType, T, MEM>::operator()(
        AMRStateData& a_output,
        const AMRStateData& a_state,
        const AMRAuxData&   a_aux,
        T                   a_scale) const
{
    // assumes boundary data is already interpolated
    // and averaged down
    
    //TODO: Check for grid equavalence
    for (int lvl = 0; lvl < m_grid.numLevels(); lvl++)
    {
        auto& output_lvl = a_output[lvl];
        auto& state_lvl  = a_state[lvl];
        auto& aux_lvl    = a_aux[lvl];
        auto& op_lvl     = m_levelOps[lvl];
        
        op_lvl(output_lvl, state_lvl, aux_lvl, a_scale);
    }

    //TODO: It's possible to reuse the flux from the apply to do refluxing
    //      with a little bit of careful reorganizing. -CLG
    for (int lvl = 1; lvl < m_grid.numLevels(); lvl++)
    {
        auto& output_crse = a_output[lvl-1];

        auto& state_crse  = a_state[lvl-1];
        auto& aux_crse    = a_aux[lvl-1];
        auto& op_crse     = m_levelOps[lvl-1];

        auto& state_fine  = a_state[lvl];
        auto& aux_fine    = a_aux[lvl];
        auto& op_fine     = m_levelOps[lvl];
        
        auto& fr_crse = *(m_fluxRegisters[lvl-1]);
        
        fr_crse.reset();
        op_crse.coarseReflux(fr_crse, state_crse, aux_crse);
        op_fine.fineReflux(fr_crse, state_fine, aux_fine);
        
        double scale = a_scale*op_crse.fluxScale() / m_dx[lvl-1];
        fr_crse.reflux(output_crse, scale);
    }
}

template <template<typename, MemType> class OPType, typename T, MemType MEM>
void AMROp<OPType, T, MEM>::operator()(
        AMRStateData&       a_output,
        const AMRStateData& a_state,
        T                   a_scale) const
{
    // assumes boundary data is already interpolated
    // and averaged down
    
    //TODO: Check for grid equavalence
    for (int lvl = 0; lvl < m_grid.numLevels(); lvl++)
    {
        auto& output_lvl = a_output[lvl];
        auto& state_lvl  = a_state[lvl];
        auto& op_lvl     = m_levelOps[lvl];
        
        op_lvl(output_lvl, state_lvl, a_scale);
    }

    //TODO: It's possible to reuse the flux from the apply to do refluxing
    //      with a little bit of careful reorganizing. -CLG
    for (int lvl = 1; lvl < m_grid.numLevels(); lvl++)
    {
        auto& output_crse = a_output[lvl-1];

        auto& state_crse  = a_state[lvl-1];
        auto& op_crse     = m_levelOps[lvl-1];

        auto& state_fine  = a_state[lvl];
        auto& op_fine     = m_levelOps[lvl];
        
        auto& fr_crse = *(m_fluxRegisters[lvl-1]);
        
        fr_crse.reset();
        op_crse.coarseReflux(fr_crse, state_crse);
        op_fine.fineReflux(fr_crse, state_fine);
        
        double scale = a_scale*op_crse.fluxScale() / m_dx[lvl-1];
        fr_crse.reflux(output_crse, scale);
    }
}

template <template<typename, MemType> class OPType, typename T, MemType MEM>
void AMROp<OPType, T, MEM>::levelApply(
        LevelStateData& a_output,
        AMRStateData& a_state, 
        AMRAuxData& a_aux, 
        int a_level,
        T   a_scale) const
{
    // assumes boundary data is already interpolated
    // and averaged down
    
    //TODO: Check for grid equavalence

    auto& state = a_state[a_level];
    auto& aux   = a_aux[a_level];
    auto& lvlOp = m_levelOps[a_level];
     
    lvlOp(a_output, state, aux, a_scale);

    if (a_level < a_state.numLevels() - 1)
    {
        auto& fineState = a_state[a_level+1];
        auto& fineAux   = a_aux[a_level+1];
        auto& fluxRegister = *(m_fluxRegisters[a_level]);
        auto& fineOp = m_levelOps[a_level+1];
         
        fluxRegister.reset();
        lvlOp.coarseReflux(fluxRegister, state, aux);
        fineOp.fineReflux(fluxRegister, fineState, fineAux);

        double scale = a_scale*lvlOp.fluxScale() / m_dx[a_level];
        fluxRegister.reflux(a_output, scale);
    }
}


template <template<typename, MemType> class OPType, typename T, MemType MEM>
void AMROp<OPType, T, MEM>::levelApply(
        LevelStateData& a_output,
        AMRStateData& a_state, 
        int a_level,
        T   a_scale) const
{
    // assumes boundary data is already interpolated
    // and averaged down
    
    //TODO: Check for grid equavalence
    
    auto& crseState = a_state[a_level];
    auto& crseOp = m_levelOps[a_level];

    crseState.exchange();    
    crseOp(a_output, crseState, a_scale);
    if (a_level < a_state.numLevels() - 1)
    {
        auto& fineState = a_state[a_level+1];
        fineState.exchange();
        reflux(a_output, crseState, fineState, a_level, a_scale);
    }
}

template <template<typename, MemType> class OPType, typename T, MemType MEM>
void AMROp<OPType, T, MEM>::reflux(
        LevelStateData& a_crseOut,
        LevelStateData& a_crseState,
        LevelStateData& a_fineState,
        int a_crseLevel,
        T   a_scale) const
{
    //std::cout << "AMROp::reflux" << std::endl;
    auto& fluxRegister = *(m_fluxRegisters[a_crseLevel]);
    auto& crseOp = m_levelOps[a_crseLevel];
    auto& fineOp = m_levelOps[a_crseLevel+1];
    fluxRegister.reset();
    crseOp.coarseReflux(fluxRegister, a_crseState);
    fineOp.fineReflux(fluxRegister, a_fineState);
    double scale = a_scale * crseOp.fluxScale() / m_dx[a_crseLevel];
    LevelStateData tmp(a_crseOut.layout(), Point::Zeros());
    tmp.setToZero();
    fluxRegister.reflux(tmp, scale);
    //std::cout << "Integral of reflux correction: " << tmp.integrate(m_dx[a_crseLevel]) << std::endl;
    //std::cout << "Integral of data BEFORE reflux: " << a_crseOut.integrate(m_dx[a_crseLevel]) << std::endl;
    fluxRegister.reflux(a_crseOut, scale);
    //std::cout << "Integral of data AFTER reflux: " << a_crseOut.integrate(m_dx[a_crseLevel]) << std::endl;
}   


template <template<typename, MemType> class OPType, typename T, MemType MEM>
void AMROp<OPType, T, MEM>::setDiagScale(double a_value, int a_gridFactor)
{
    double value = a_value;
    double factor = pow(PR_AMR_REFRATIO, a_gridFactor);
    for (int lvl = 0; lvl < m_grid.numLevels(); lvl++)
    {
        m_levelOps[lvl].setDiagScale(value);
        value /= factor;
    }
}

template <template<typename, MemType> class OPType, typename T, MemType MEM>
void AMROp<OPType, T, MEM>::setFluxScale(double a_value, int a_gridFactor)
{
    double value = a_value;
    double factor = pow(PR_AMR_REFRATIO, a_gridFactor);
    for (int lvl = 0; lvl < m_grid.numLevels(); lvl++)
    {
        m_levelOps[lvl].setFluxScale(value);
        value /= factor;
    }
}
