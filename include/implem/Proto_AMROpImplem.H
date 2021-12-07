
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
        const AMRAuxData&   a_aux) const
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
        
        op_lvl(output_lvl, state_lvl, aux_lvl);
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
        
        op_crse.coarseReflux(fr_crse, state_crse, aux_crse);
        op_fine.fineReflux(fr_crse, state_fine, aux_fine);
        
        double scale = op_crse.fluxScale() / m_dx[lvl-1];
        fr_crse.reflux(output_crse, scale);
    }
}

template <template<typename, MemType> class OPType, typename T, MemType MEM>
void AMROp<OPType, T, MEM>::operator()(
        AMRStateData& a_output,
        const AMRStateData& a_state) const
{
    // assumes boundary data is already interpolated
    // and averaged down
    
    //TODO: Check for grid equavalence
    for (int lvl = 0; lvl < m_grid.numLevels(); lvl++)
    {
        auto& output_lvl = a_output[lvl];
        auto& state_lvl  = a_state[lvl];
        auto& op_lvl     = m_levelOps[lvl];
        
        op_lvl(output_lvl, state_lvl);
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
        
        op_crse.coarseReflux(fr_crse, state_crse);
        op_fine.fineReflux(fr_crse, state_fine);
        
        double scale = op_crse.fluxScale() / m_dx[lvl-1];
        fr_crse.reflux(output_crse, scale);
    }
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
