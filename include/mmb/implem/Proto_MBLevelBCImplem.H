template<typename T, unsigned int C, 
    template<MemType> typename MAP,
    MemType MEM, Centering CTR>
MBLevelBC<T,C,MAP,MEM,CTR>::MBLevelBC(const MBLevelMap<MAP,MEM>& a_map)
{
    m_map = a_map;
    this->init();
}
