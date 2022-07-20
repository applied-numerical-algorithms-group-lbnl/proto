
template<typename T, unsigned int C, MemType MEM, Centering CTR>
MBLevelBoxData<T, C, MEM, CTR>::MBLevelBoxData(
        MBDisjointBoxLayout& a_layout,
        unsigned int a_ghost)
{
    std::array<Point, DIM+1> ghost;
    ghost.fill(Point::Ones(a_ghost));
    define(a_layout, ghost);
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
MBLevelBoxData<T, C, MEM, CTR>::MBLevelBoxData(
        MBDisjointBoxLayout&      a_layout,
        std::array<Point, DIM+1> a_ghost)
{
    define(a_layout, a_ghost);
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void MBLevelBoxData<T, C, MEM, CTR>::define(
        MBDisjointBoxLayout&      a_layout,
        std::array<Point, DIM+1> a_ghost)
{
    // a_ghost convention:
    // [interior, codim1, codim2, ...]
    m_ghost = a_ghost;
    m_layout = &a_layout;
    m_data.clear();
    m_data.resize(a_layout.numBlocks());
    for (int bi = 0; bi < a_layout.numBlocks(); bi++)
    {
        m_data[bi] = std::make_shared<LevelBoxData<T, C, MEM, CTR>>(
                a_layout.layout(bi), a_ghost[0]);
    }
    for (auto iter : a_layout)
    {
        auto& patch = (*this)[iter];
    }
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
template<typename Func, typename... Srcs>
void MBLevelBoxData<T, C, MEM, CTR>::initialize(Func& a_func, Srcs... a_srcs);

template<typename T, unsigned int C, MemType MEM, Centering CTR>
template<typename Func, typename... Srcs>
void MBLevelBoxData<T, C, MEM, CTR>::initConvolve(Func& a_func, Srcs... a_srcs);

template<typename T, unsigned int C, MemType MEM, Centering CTR>
BoxData<T, C, MEM>& MBLevelBoxData<T, C, MEM, CTR>::operator[](index_t a_index)
{
    unsigned int block = m_layout.block(a_index);
    auto localIndex = m_layout.blockIndex[a_index];
    return m_data[block][localIndex];
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
const BoxData<T, C, MEM>& MBLevelBoxData<T, C, MEM, CTR>::operator[](index_t a_index) const
{
    unsigned int block = m_layout.block(a_index);
    auto localIndex = m_layout.blockIndex[a_index];
    return m_data[block][localIndex];
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
LevelBoxData<T, C, MEM, CTR>& MBLevelBoxData<T, C, MEM, CTR>::block(unsigned int a_block)
{
    return m_data[a_block];
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
const LevelBoxData<T, C, MEM, CTR>&
MBLevelBoxData<T, C, MEM, CTR>::block(unsigned int a_block) const
{
    return m_data[a_block];
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
iterator_t MBLevelBoxData<T, C, MEM, CTR>::begin()
{
    return m_layout->begin();
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
iterator_t MBLevelBoxData<T, C, MEM, CTR>::end()
{
    return m_layout->end();
}
