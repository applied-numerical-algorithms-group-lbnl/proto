template<typename T, unsigned int C,
    typename MAP,
    MemType MEM, Centering CTR>
void ConstBC<T,C,MAP,MEM,CTR>::setVal(T a_value, int a_comp)
{
    if (a_comp < 0)
    {
        m_compValues.clear();
        m_globalValue = a_value;
    } else {
        m_compValues.push_back(std::make_pair(a_comp, a_value));
    }
}

template<typename T, unsigned int C,
    typename MAP,
    MemType MEM, Centering CTR>
void ConstBC<T,C,MAP,MEM,CTR>::apply(MBLevelBoxData<T,C,MEM,CTR>& a_data) const
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
