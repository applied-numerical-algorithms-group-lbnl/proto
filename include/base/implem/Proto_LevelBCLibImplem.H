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
}

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



