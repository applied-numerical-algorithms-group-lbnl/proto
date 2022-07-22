
template<typename T, unsigned int C, MemType MEM, Centering CTR>
MBLevelBoxData<T, C, MEM, CTR>::MBLevelBoxData(
        MBDisjointBoxLayout& a_layout,
        unsigned int a_ghost)
{
    std::array<unsigned int, DIM+1> ghost;
    ghost.fill(a_ghost);
    define(a_layout, ghost);
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
MBLevelBoxData<T, C, MEM, CTR>::MBLevelBoxData(
        MBDisjointBoxLayout&      a_layout,
        std::array<unsigned int, DIM+1> a_ghost)
{
    define(a_layout, a_ghost);
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void MBLevelBoxData<T, C, MEM, CTR>::define(
        MBDisjointBoxLayout&      a_layout,
        std::array<unsigned int, DIM+1> a_ghost)
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
                a_layout.layout(bi), Point::Ones(a_ghost[0]));
    }
    for (int bi = 0; bi < a_layout.numBlocks(); bi++)
    {
        Box K = Box::Kernel(1);
        Box blockDomainBox = a_layout.domain().blockDomain(bi).box();
        for (auto di : K)
        {
            if (di == Point::Zeros()) { continue; }
            unsigned int boundCodim = codim(di);
            boundsKey_t key(bi, di);
            Box localBoundBox = blockDomainBox.adjacent(di, a_ghost[boundCodim]);
            auto bounds = a_layout.domain().graph().boundaries(bi, di);
            for (auto bound : bounds)
            {
                auto adjBlock = bound.dstBlock;
                Point diAdj = bound.dstToSrc;
                Box adjBlockDomainBox = a_layout.domain().blockDomain(bound.dstBlock).box();
                Box adjBoundBox = adjBlockDomainBox.edge(diAdj, a_ghost[boundCodim]);
                
                MBBoundaryData<T, C, MEM> boundData(
                        bi, adjBlock, localBoundBox, adjBoundBox);
                m_bounds[key].push_back(boundData);
            }
        }
    }
/*
    for (auto iter : a_layout)
    {
        auto block = a_layout.block(iter);
        auto blockLayout = a_layout.layout(block);
        auto blockIndex = a_layout.blockIndex(iter);
        auto& patch = (*this)[iter];
        Box patchBox = blockLayout[blockIndex];
        Point patchID = blockLayout.point(blockIndex);
        Box blockPatches = blockLayout.patchDomain().box();
        Box K = Box::Kernel(1);
        for (auto di : K)
        {
            if (di == Point::Zeros()) { continue; }
            Box blockDomainBoundary = blockPatches.edge(di, 1);
            if (blockDomainBoundary.contains(patchID))
            {
                unsigned int boundCodim = codim(di);
                Box localBoundBox = patchBox.adjacent(di,a_ghost[boundCodim]);
                auto bounds = a_layout.domain().graph().boundaries(block, di);
                for (auto bound : bounds)
                {
                    auto adjBlock = bound.dstBlock;
                    Point diAdj = bound.dstToSrc;
                    Box adjBlockDomainBox = a_layout.domain().blockDomain(bound.dstBlock).box();
                    Box adjBoundBox = adjBlockDomainBox.edge(diAdj, a_ghost[boundCodim]);
                    boundsKey_t key(di, iter);
                    MBBoundaryData<T, C, MEM> boundData(
                            block, adjBlock, localBoundBox, adjBoundBox);
                    m_bounds[key].push_back(boundData);
                }
            }
        }
    }
    */
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
template<typename Func, typename... Srcs>
void MBLevelBoxData<T, C, MEM, CTR>::initialize(Func& a_func, Srcs... a_srcs) {}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
template<typename Func, typename... Srcs>
void MBLevelBoxData<T, C, MEM, CTR>::initConvolve(Func& a_func, Srcs... a_srcs){}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
BoxData<T, C, MEM>&
MBLevelBoxData<T, C, MEM, CTR>::operator[](index_t a_index)
{
    unsigned int block = m_layout->block(a_index);
    auto localIndex = m_layout->blockIndex(a_index);
    return (*m_data[block])[localIndex];
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
const BoxData<T, C, MEM>&
MBLevelBoxData<T, C, MEM, CTR>::operator[](index_t a_index) const
{
    unsigned int block = m_layout->block(a_index);
    auto localIndex = m_layout->blockIndex(a_index);
    return (*m_data[block])[localIndex];
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
LevelBoxData<T, C, MEM, CTR>&
MBLevelBoxData<T, C, MEM, CTR>::blockData(unsigned int a_block)
{
    return m_data[a_block];
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
const LevelBoxData<T, C, MEM, CTR>&
MBLevelBoxData<T, C, MEM, CTR>::blockData(unsigned int a_block) const
{
    return m_data[a_block];
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
std::vector<MBBoundaryData<T, C, MEM>>&
MBLevelBoxData<T, C, MEM, CTR>::bounds(unsigned int a_block, Point a_dir)
{
    return m_bounds[boundsKey_t(a_block, a_dir)];
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
const std::vector<MBBoundaryData<T, C, MEM>>&
MBLevelBoxData<T, C, MEM, CTR>::bounds(unsigned int a_block, Point a_dir) const
{
    return m_bounds[boundsKey_t(a_block, a_dir)];
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
typename MBLevelBoxData<T, C, MEM, CTR>::iterator_t
MBLevelBoxData<T, C, MEM, CTR>::begin()
{
    return m_layout->begin();
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
typename MBLevelBoxData<T, C, MEM, CTR>::iterator_t
MBLevelBoxData<T, C, MEM, CTR>::end()
{
    return m_layout->end();
}
