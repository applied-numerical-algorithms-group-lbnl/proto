
template<typename T, unsigned int C, MemType MEM, Centering CTR>
void AMRData<T, C, MEM, CTR>::define(AMRGrid& a_grid, Point a_ghost)
{
    m_ghost = a_ghost;
    m_grid = a_grid;
    m_data.clear();
    for (int ii = 0; ii < a_grid.numLevels(); ii++)
    {
        auto level = std::make_shared<LevelBoxData<T, C, MEM, CTR>>(a_grid[ii], a_ghost);
        m_data.push_back(level);
    }

    // interpolation is 5th order by default. 
    // TODO: This order should be controlled somewhere else (compile time constant?)
    m_boundInterpOrder = -1;
    buildInterpStencils(5);
    m_defined = true;
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
LevelBoxData<T, C, MEM, CTR>&
AMRData<T, C, MEM, CTR>::operator[](unsigned int a_level)
{
    PROTO_ASSERT(a_level < m_data.size(),
        "AMRData::operator[] | Error: level %u is out of bounds.", a_level);
    return *m_data[a_level];
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
const LevelBoxData<T, C, MEM, CTR>&
AMRData<T, C, MEM, CTR>::operator[](unsigned int a_level) const
{
    PROTO_ASSERT(a_level < m_data.size(),
        "AMRData::operator[] | Error: level %u is out of bounds.", a_level);
    return *m_data[a_level];
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void
AMRData<T, C, MEM, CTR>::setToZero()
{
    for (int ii = 0; ii < m_data.size(); ii++)
    {
        m_data[ii]->setToZero();
    }
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
template<typename Func, typename... Srcs>
void
AMRData<T, C, MEM, CTR>::initialize(double a_dx0, Func& a_func, Srcs... a_srcs)
{
    double dx = a_dx0;
    for (int ii = 0; ii < m_data.size(); ii++)
    {
        m_data[ii]->initialize(a_func, dx, a_srcs...);
        dx /= PR_AMR_REFRATIO;
    }
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
template<typename Func, typename... Srcs>
void
AMRData<T, C, MEM, CTR>::initConvolve(double a_dx0, Func& a_func, Srcs... a_srcs)
{
    double dx = a_dx0;
    for (int ii = 0; ii < m_data.size(); ii++)
    {
        m_data[ii]->initConvolve(a_func, dx, a_srcs...);
        dx /= PR_AMR_REFRATIO;
    }
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void
AMRData<T, C, MEM, CTR>::averageDown()
{
    if (numLevels() < 2) { return; }
    auto average = Stencil<double>::AvgDown(PR_AMR_REFRATIO);
    for (int lvl = numLevels() - 2; lvl >= 0; lvl--)
    {
        auto& crse = operator[](lvl);
        auto& fine = operator[](lvl+1);
        auto  cfLayout = fine.layout().coarsen(Point::Ones(PR_AMR_REFRATIO));
        LevelBoxData<T, C, MEM, CTR> temp(cfLayout, Point::Zeros());
        for (auto iter = fine.layout().begin(); iter.ok(); ++iter)
        {
            auto& fine_i = fine[*iter];
            auto& crse_i = temp[*iter];
            
            crse_i |= average(fine_i);
        }
        temp.copyTo(crse);
    }
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void
AMRData<T, C, MEM, CTR>::interp(int a_order)
{
    //FIXME: This implementation doesn't address refined regions directly
    //       adjacent to domain boundaries. -CLG
    if (numLevels() < 2) { return; }
    buildInterpStencils(a_order);
    for (int lvl = 1; lvl < numLevels(); lvl++)
    {
        auto& crse = operator[](lvl-1);
        auto& fine = operator[](lvl);
        interpBoundaries(crse, fine, m_boundInterp);
    }
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void
AMRData<T, C, MEM, CTR>::interpLevel(int a_level, int a_order)
{
    //FIXME: This implementation doesn't address refined regions directly
    //       adjacent to domain boundaries. -CLG
    if (numLevels() < 2) { return; }
    buildInterpStencils(a_order);

    auto& crse = operator[](a_level-1);
    auto& fine = operator[](a_level);
    interpBoundaries(crse, fine, m_boundInterp);
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void
AMRData<T, C, MEM, CTR>::buildInterpStencils(int a_order)
{
    if (a_order == m_boundInterpOrder) { return; }
    m_boundInterpOrder = a_order;
    int shiftMax;
    Box kernel;
    int polyOrder;
    switch (a_order)
    {
        case 3:
        {
            shiftMax = 2;
            kernel = Box::Kernel(2);
            polyOrder = 2;
            break;
        }
        case 5:
        {
            shiftMax = 4;
            kernel = Box::Kernel(2);
            polyOrder = 4;
            break;
        }
        default:
        {
            MayDay<void>::Error(
                "AMRData::buildInterpStencils | Error: Invalid order. Valid orders: 3, 5.");
            break;
        }
    }
    m_boundInterp = InterpStencil<double>::Build(shiftMax, kernel, polyOrder, PR_AMR_REFRATIO);
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
double
AMRData<T, C, MEM, CTR>::integrate(double a_cdx, unsigned int a_c) const
{
    double integral = operator[](0).integrate(a_cdx, a_c);
   
    // technically, so long as averageDown() was called before integrate,
    // the rest of this code does nothing.
    /* 
    double cdx = a_cdx;
    double fdx = cdx / PR_AMR_REFRATIO;

    for (int lvl = 1; lvl < numLevels(); lvl++)
    {
        auto& crse = operator[](lvl-1);
        auto& fine = operator[](lvl);

        auto cfLayout = m_grid[lvl].coarsen(Point::Ones(PR_AMR_REFRATIO));
        LevelBoxData<T, C, MEM, CTR> cfTemp(cfLayout, Point::Zeros());
        //TODO: improve performance using componentwise copy (we only need a_c)
        crse.copyTo(cfTemp);
        integral -= cfTemp.integrate(cdx, a_c);
        integral += fine.integrate(fdx, a_c);

        cdx /= PR_AMR_REFRATIO;
        fdx /= PR_AMR_REFRATIO;
    }
    */
    return integral;
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
double
AMRData<T, C, MEM, CTR>::absMax(unsigned int a_c) const
{
    return operator[](0).absMax(a_c);
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void
AMRData<T, C, MEM, CTR>::exchange()
{
    for (int lvl = 0; lvl < numLevels(); lvl++)
    {
        m_data[lvl]->exchange();
    }
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void
AMRData<T, C, MEM, CTR>::regrid(unsigned int a_level)
{
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void
AMRData<T, C, MEM, CTR>::computeTags(
    LevelTagData& a_tags,
    unsigned int a_level,
    Point a_bufferSize,
    T a_threshold) const
{
    PROTO_ASSERT(a_level != maxLevels()-1,
        "AMRData::computeTags | Error: Cannot generate tags on the finest level");
    PROTO_ASSERT(a_level < maxLevels()-1,
        "AMRData::computeTags | Error: Level %u is out of bounds", a_level);
    PROTO_ASSERT(a_level < m_data.size(),
        "AMRData::computeTags | Error: No data defined on level %u.", a_level);

    a_tags.define(m_grid[a_level], a_bufferSize);

    std::vector<Stencil<T>> D;
    for (int dir = 0; dir < DIM; dir++)
    {
        Stencil<T> Si = ((T)1)*Shift::Basis(dir, 1) - ((T)1)*Shift::Basis(dir, -1);
        D.push_back(Si);
    }

    for (auto iter = m_grid[a_level].begin(); iter.ok(); ++iter)
    {
        const auto& data = operator[](a_level)[*iter];
        auto& tags = a_tags[*iter];
        BoxData<T, DIM> diffs(iter.box()); 
        for (int dir = 0; dir < DIM; dir++)
        {
            BoxData<T, 1> diff = slice(diffs, dir);
            diff |= D[dir](data);
        }
        forallInPlace(
        [] PROTO_LAMBDA (Var<char, 1>& v_tags, Var<T, DIM>& v_diffs, T v_threshold)
        {
            T diffSq = 0;
            for (int dir = 0; dir < DIM; dir++)
            {
                diffSq += v_diffs(dir)*v_diffs(dir);
            }
            T diff = sqrt(diffSq);
            if (diff > v_threshold)
            {
                v_tags(0) = 1;
            } else {
                v_tags(0) = 0;
            }
        }, tags, diffs, a_threshold);
    }
    
    AMRGrid::buffer(a_tags, a_bufferSize);
}

