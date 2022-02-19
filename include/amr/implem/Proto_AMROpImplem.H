
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
    m_interp.resize(numLevels - 1); 
    m_grid = a_grid; //intentional deep copy
    double dx = a_cdx;
    for (int lvl = 0; lvl < numLevels; lvl++)
    {
        m_dx[lvl] = dx;
        m_levelOps[lvl].define(dx);
        if (lvl < numLevels - 1)
        {
            Point refRatio = m_grid.refRatio(lvl);
            m_fluxRegisters[lvl] =
                std::make_shared<LevelFluxRegister<T>>(
                    a_grid[lvl], a_grid[lvl+1],
                    refRatio);
            m_interp[lvl] = 
                InterpStencil<T>::Build(BOP::order() + 1, refRatio);
            //FIXME: Assumes isotropic refinement
            dx /= refRatio[0];
        }
    }
    m_defined = true;
}

template <template<typename, MemType> class OPType, typename T, MemType MEM>
void AMROp<OPType, T, MEM>::operator()(
        AMRStateData& a_output,
        const AMRStateData& a_state,
        const AMRAuxData&   a_aux,
        T                   a_scale)
{
    PROTO_ASSERT(a_output.layout().compatible(a_state.layout()),
        "AMROp::operator() | Error: output and state have incompatible layouts.");
    PROTO_ASSERT(a_output.layout().compatible(a_aux.layout()),
        "AMROp::operator() | Error: output and aux have incompatible layouts.");
    
    a_state.averageDown();
    for (int lvl = 0; lvl < m_grid.numLevels(); lvl++)
    {
        auto& output = a_output[lvl];
        auto& state  = a_state[lvl];
        auto& aux    = a_aux[lvl];
        auto& op     = m_levelOps[lvl];
        
        if (lvl > 0)
        {
            auto& crseState = a_state[lvl-1];
            interpBoundaries(crseState, state, m_interp);
        }
        op(output, state, aux, a_scale);
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
        AMRStateData& a_output,
        AMRStateData& a_state,
        T             a_scale)
{
    PROTO_ASSERT(a_output.grid().compatible(a_state.grid()),
        "AMROp::operator() | Error: output and state have incompatible grids.");
    
    // average down all levels recursively
    a_state.averageDown(); 
  
    for (int lvl = 0; lvl < m_grid.numLevels(); lvl++)
    {
        auto& output = a_output[lvl];
        auto& state  = a_state[lvl];
        auto& op     = m_levelOps[lvl];
        
        // if a coarser level exists, interpolate boundaries
        if (lvl > 0)
        {
            auto& crseState = a_state[lvl-1];
            interpBoundaries(crseState, state, m_interp[lvl-1]);
        }
        
        // apply the level operator
        op(output, state, a_scale);
    }

    // compute refluxing correction
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
        
        // NB: a_scale multiplies the whole operator
        //     fluxScale multiplies only the div(flux) and source terms
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
        T   a_scale)
{
    PROTO_ASSERT(a_output.layout().compatible(a_state[a_level].layout()),
        "AMROp::operator() | Error: output and state[%i] have incompatible layouts.", a_level);
    PROTO_ASSERT(a_output.layout().compatible(a_aux[a_level].layout()),
        "AMROp::operator() | Error: output and aux[%i] have incompatible layouts.", a_level);

    // assumes boundary data is already interpolated
    // and averaged down

    auto& state = a_state[a_level];
    auto& aux   = a_aux[a_level];
    auto& op = m_levelOps[a_level];
    
    // if there is a finer level, average down to this level
    if (a_level < m_grid.numLevels() - 1)
    {
        // aux data shouldn't need averaging
        Point refRatio = m_grid.refRatio(a_level);
        auto& fineState = a_state[a_level + 1];
        averageDown(state, fineState, refRatio);
    }
    
    // if there is a coarser level, average down to it and then
    // interpolate boundary conditions from it
    if (a_level > 0)
    {
        // aux data shouldn't need averaging / interpolation
        Point refRatio = m_grid.refRatio(a_level-1);
        auto& crseState = a_state[a_level-1];
        averageDown(crseState, state, refRatio);
        interpBoundaries(crseState, state, m_interp[a_level-1]);
    }
     
    op(a_output, state, aux, a_scale);

    // TODO: write this as a call to AMROp::reflux
    if (a_level < a_state.numLevels() - 1)
    {
        auto& fineState = a_state[a_level+1];
        auto& fineAux   = a_aux[a_level+1];
        auto& fluxRegister = *(m_fluxRegisters[a_level]);
        auto& fineOp = m_levelOps[a_level+1];
         
        fluxRegister.reset();
        op.coarseReflux(fluxRegister, state, aux);
        fineOp.fineReflux(fluxRegister, fineState, fineAux);

        double scale = a_scale*op.fluxScale() / m_dx[a_level];
        fluxRegister.reflux(a_output, scale);
    }
}


template <template<typename, MemType> class OPType, typename T, MemType MEM>
void AMROp<OPType, T, MEM>::levelApply(
        LevelStateData& a_output,
        AMRStateData& a_state, 
        int a_level,
        T   a_scale)
{
    PROTO_ASSERT(a_output.layout().compatible(a_state[a_level].layout()),
        "AMROp::operator() | Error: output and state[%i] have incompatible layouts.", a_level);
     
    auto& state = a_state[a_level];
    auto& op = m_levelOps[a_level];
    
    // if there is a finer level, average down to this level
    if (a_level < m_grid.numLevels() - 1)
    {
        Point refRatio = m_grid.refRatio(a_level);
        auto& fineState = a_state[a_level + 1];
        averageDown(state, fineState, refRatio);
    }
    
    // if there is a coarser level, average down to it and then
    // interpolate boundary conditions from it
    if (a_level > 0)
    {
        Point refRatio = m_grid.refRatio(a_level-1);
        auto& crseState = a_state[a_level-1];
        averageDown(crseState, state, refRatio);
        interpBoundaries(crseState, state, m_interp[a_level-1]);
    }
    
    // apply the level operator
    op(a_output, state, a_scale);

    // if there is a finer level, use it to compute reflux correction
    if (a_level < a_state.numLevels() - 1)
    {
        auto& fineState = a_state[a_level + 1];
        reflux(a_output, state, fineState, a_level, a_scale);
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
    auto& fluxRegister = *(m_fluxRegisters[a_crseLevel]);
    auto& crseOp = m_levelOps[a_crseLevel];
    auto& fineOp = m_levelOps[a_crseLevel+1];
    fluxRegister.reset();
    crseOp.coarseReflux(fluxRegister, a_crseState);
    fineOp.fineReflux(fluxRegister, a_fineState);
    double scale = a_scale * crseOp.fluxScale() / m_dx[a_crseLevel];
    fluxRegister.reflux(a_crseOut, scale);
}   

template <template<typename, MemType> class OPType, typename T, MemType MEM>
InterpStencil<T>& AMROp<OPType, T, MEM>::interpStencil(int a_level)
{
    PROTO_ASSERT(a_level < m_grid.numLevels()-1,
        "AMROp::interpStencil | Error: level %i is out of bounds.", a_level);
    return m_interp[a_level];
}

template <template<typename, MemType> class OPType, typename T, MemType MEM>
void AMROp<OPType, T, MEM>::setDiagScale(double a_value, int a_gridFactor)
{
    double value = a_value;
    for (int lvl = 0; lvl < m_grid.numLevels(); lvl++)
    {
        m_levelOps[lvl].setDiagScale(value);
        //FIXME: Assumes isotropic refinement
        if (lvl < m_grid.numLevels() - 1)
        {
            Point refRatio = m_grid.refRatio(lvl);
            double factor = pow(refRatio[0], a_gridFactor);
            value /= factor;
        }
    }
}

template <template<typename, MemType> class OPType, typename T, MemType MEM>
void AMROp<OPType, T, MEM>::setFluxScale(double a_value, int a_gridFactor)
{
    double value = a_value;
    for (int lvl = 0; lvl < m_grid.numLevels(); lvl++)
    {
        m_levelOps[lvl].setFluxScale(value);
        //FIXME: Assumes isotropic refinement
        if (lvl < m_grid.numLevels() - 1)
        {
            Point refRatio = m_grid.refRatio(lvl);
            double factor = pow(refRatio[0], a_gridFactor);
            value /= factor;
        }
    }
}
