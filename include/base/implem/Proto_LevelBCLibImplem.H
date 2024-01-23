template<typename T, unsigned int C, MemType MEM, Centering CTR>
void PeriodicBC<T,C,MEM,CTR>::apply(LevelBoxData<T,C,MEM,CTR>& a_data) const
{
    a_data.exchange();
}

//TODO: This isn't quite right; it sets ghost values instead of doing extrapolation
template<typename T, unsigned int C, MemType MEM, Centering CTR>
void ConstDirichletBC<T,C,MEM,CTR>::apply(LevelBoxData<T,C,MEM,CTR>& a_data) const
{
    if (m_compValues.size() == 0)
    {
        a_data.setBoundary(m_globalValue);
    } else {
        for (auto val : m_compValues)
        {
            a_data.setBoundary(val.second, val.first);
        }
    }
    a_data.exchange(); // overwrite the periodic cells
    
    // this works, but it's probably slower
#if 0
    for (auto iter : a_data.layout())
    {
        for (auto dir : Box::Kernel(1))
        {
            std::cout << "patch: " << a_data.layout().point(iter);
            std::cout << " | dir: " << dir << std::endl;
            if (a_data.layout().isDomainBoundary(iter, dir))
            {
                Box boundBox = a_data.layout()[iter].adjacent(dir*a_data.ghost());
                std::cout << "\tfound boundary with box: " << boundBox << std::endl;
                auto& patch = a_data[iter];
                PROTO_ASSERT(patch.box().contains(boundBox),
                        "ConstDirichletBC::apply | Error: data does not contain boundary");
                forallInPlace(
                        [] PROTO_LAMBDA (Var<double, C, MEM>& v_data, double v_value)
                        {
                            for (int cc = 0; cc < C; cc++)
                            {
                                v_data(cc) = v_value;
                            }
                        }, patch, m_value);
            }
        }
    }
#endif
}

#ifdef PR_MMB
template<typename T, unsigned int C, MemType MEM, Centering CTR>
void ConstDirichletBC<T,C,MEM,CTR>::apply(MBLevelBoxData<T,C,MEM,CTR>& a_data) const
{
    if (m_compValues.size() == 0)
    {
        a_data.setBoundary(m_globalValue);
    } else {
        for (auto val : m_compValues)
        {
            a_data.setBoundary(val.second, val.first);
        }
    }
}
#endif
template<typename T, unsigned int C, MemType MEM, Centering CTR>
void ConstDirichletBC<T,C,MEM,CTR>::setVal(T a_value, int a_comp)
{
    if (a_comp < 0)
    {
        m_compValues.clear();
        m_globalValue = a_value;
    } else {
        m_compValues.push_back(std::make_pair(a_comp, a_value));
    }
}



