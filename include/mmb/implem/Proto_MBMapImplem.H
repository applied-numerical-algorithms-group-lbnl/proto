#if 1

template<typename Func, MemType MEM>
MBMap<Func, MEM>::MBMap(
        const Func& a_func,
        MBDisjointBoxLayout& a_layout,
        Array<Point, DIM+1>& a_ghost,
        unsigned int a_order)
{
    define(a_func, a_layout, a_ghost, a_order);
}

template<typename Func, MemType MEM>
void MBMap<Func, MEM>::define(
        const Func& a_func,
        MBDisjointBoxLayout& a_layout,
        Array<Point, DIM+1>& a_ghost,
        unsigned int a_order)
{
    m_func = &a_func;
    m_order = a_order;
    m_X.define(a_layout, a_ghost);
    for (auto iter : a_layout)
    {
        auto& patch = m_X[iter];
        auto block = a_layout.block(iter);
        Point blockSize = a_layout.blockLayout(block).domain().box().sizes();
        forallInPlace_p(*m_func, patch, block, blockSize);
    }
}

template<typename Func, MemType MEM>
MBLevelBoxData<double, 3, MEM, PR_NODE>&
MBMap<Func, MEM>::map()
{
    return m_X;
}

template<typename Func, MemType MEM>
bool MBMap<Func, MEM>::contains(
        const MBIndex& a_patch,
        const Point& a_cell,
        const Array<double, 3>& a_X) const
{
    std::cout << "MBMap::contains is not yet implemented!!" << std::endl;
    //const auto& layout = m_X.layout();
    //auto const& X = m_X[a_patch];
    //Box b = X.box();
    //std::array<T, 3> x0, x1;
    return true;
}

#else
    template<typename T, MemType MEM>
MBMap<T, MEM>::MBMap(
        MBDisjointBoxLayout& a_layout,
        std::array<Point, DIM+1>& a_ghost,
        unsigned int a_order)
{
    //TODO: will need an extra layer of ghost cells for taking derivatives
    define(a_layout, a_ghost, a_order);
}

template<typename T, MemType MEM>
void
MBMap<T, MEM>::define(
        MBDisjointBoxLayout& a_layout,
        std::array<Point, DIM+1>& a_ghost,
        unsigned int a_order)
{
    std::array<Point, DIM+1> zeros;
    zeros.fill(Point::Zeros());
    m_X.define(a_layout, a_ghost);
    m_J.define(a_layout, zeros);
    m_order = a_order;
    m_defined = true;
    m_computed = false;
}

template<typename T, MemType MEM>
template<typename Func, typename... Srcs> 
void
MBMap<T, MEM>::compute(const Func& a_func, Srcs&&... a_srcs)
{
    PROTO_ASSERT(m_defined,
            "MBMap::compute | Error: MBMap data holders not defined.");
    for (auto iter : m_X.layout())
    {
        auto& x = m_X[iter];
        auto block = m_X.layout().block(iter);
        forallInPlace_p(a_func, x, block, a_srcs...);
        //TODO: compute Jacobian
    }
    m_X.fillBoundaries();
    m_computed = true;
}

/// Map Access
template<typename T, MemType MEM>
MBLevelBoxData<T, 3, MEM, PR_NODE>&
MBMap<T, MEM>::map()
{
    PROTO_ASSERT(m_defined,
            "MBMap::map | Error: MBMap data holders not defined.");
    PROTO_ASSERT(m_computed,
            "MBMap::map | Error: MBMap not yet computed.");
    return m_X;
}

/// Jacobian Access
template<typename T, MemType MEM>
const MBLevelBoxData<T, 1, MEM, PR_CELL>&
MBMap<T, MEM>::jacobian() const
{
    PROTO_ASSERT(m_defined,
            "MBMap::jacobian | Error: MBMap data holders not defined.");
    PROTO_ASSERT(m_computed,
            "MBMap::jacobian | Error: MBMap not yet computed.");
    return m_J;
}

/// Compute Deformation Gradient
template<typename T, MemType MEM>
BoxData<T, 3, MEM, 3>
MBMap<T, MEM>::deformation(
        unsigned int a_blockIndex,
        const LevelIndex& a_patchIndex,
        int a_dir) const
{
    PROTO_ASSERT(m_defined,
            "MBMap::deformation | Error: MBMap data holders not defined.");
    PROTO_ASSERT(m_computed,
            "MBMap::deformation | Error: MBMap not yet computed.");
}

template<typename T, MemType MEM>
BoxData<T, 3, MEM, 3>
MBMap<T, MEM>::deformation(
        const MBIndex& a_index,
        int a_dir) const
{
    PROTO_ASSERT(m_defined,
            "MBMap::deformation | Error: MBMap data holders not defined.");
    PROTO_ASSERT(m_computed,
            "MBMap::deformation | Error: MBMap not yet computed.");
}

/// Compute Cofactor Matrix
template<typename T, MemType MEM>
BoxData<T, 3, MEM, 3>
MBMap<T, MEM>::cofactors(
        unsigned int a_blockIndex,
        const LevelIndex& a_patchIndex, int a_dir) const
{
    PROTO_ASSERT(m_defined,
            "MBMap::cofactors | Error: MBMap data holders not defined.");
    PROTO_ASSERT(m_computed,
            "MBMap::cofactors | Error: MBMap not yet computed.");
}

template<typename T, MemType MEM>
BoxData<T, 3, MEM, 3>
MBMap<T, MEM>::cofactors(const MBIndex& a_index, int a_dir) const
{
    PROTO_ASSERT(m_defined,
            "MBMap::cofactors | Error: MBMap data holders not defined.");
    PROTO_ASSERT(m_computed,
            "MBMap::cofactors | Error: MBMap not yet computed.");
}


// Find Cell
template<typename T, MemType MEM>
bool
MBMap<T, MEM>::contains(
        const MBIndex& a_patch,
        const Point& a_cell,
        const std::array<T, 3>& a_X) const
{
    const auto& layout = m_X.layout();
    auto const& X = m_X[a_patch];
    Box b = X.box();
    std::array<T, 3> x0, x1;
    return true;
}
#endif // end Version 1;

/*
template<typename T, MemType MEM>
std::vector<MBIndex>
MBMap<T, MEM>::find(const std::array<T, 3>& a_X) const
{
}

template<typename T, MemType MEM>
MBIndex
MBMap<T, MEM>::find(
        const std::array<T, 3>& a_X,
        const MBIndex& a_nearPoint) const;
{
}
*/

