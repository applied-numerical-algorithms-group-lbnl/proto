
template <typename T, unsigned int C, MemType MEM, Centering CTR>
MBLevelBoxData<T, C, MEM, CTR>::MBLevelBoxData()
{
}

template <typename T, unsigned int C, MemType MEM, Centering CTR>
MBLevelBoxData<T, C, MEM, CTR>::MBLevelBoxData(
    const MBDisjointBoxLayout& a_layout,
    Point               a_ghost,
    Point               a_boundGhost)
{
    Array<Point, DIM + 1> ghost;
    ghost.fill(a_ghost);
    define(a_layout, ghost, a_boundGhost);
}

template <typename T, unsigned int C, MemType MEM, Centering CTR>
MBLevelBoxData<T, C, MEM, CTR>::MBLevelBoxData(
    const MBDisjointBoxLayout&    a_layout,
    const Array<Point, DIM+1>&    a_ghost,
    Point                   a_boundGhost)
{
    define(a_layout, a_ghost, a_boundGhost);
}

template <typename T, unsigned int C, MemType MEM, Centering CTR>
void MBLevelBoxData<T, C, MEM, CTR>::define(
    const MBDisjointBoxLayout&    a_layout,
    const Array<Point, DIM+1>&    a_ghost,
    Point                   a_boundGhost)
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
            a_layout.blockLayout(bi), a_ghost[0]);
    }
    for (auto iter : a_layout)
    {
        auto block = a_layout.block(iter);
        auto patch = a_layout.point(iter);
        auto localLayout = a_layout.blockLayout(block);
        auto patchDomain = localLayout.patchDomain().box();
        auto localIndex = localLayout.find(patch);
        auto patchBox = localLayout[localIndex];

        Box K = Box::Kernel(1);
        for (int cc = 1; cc <= DIM; cc++)
        {
            int ghostSize = m_ghost[cc].max();
            for (auto di : K)
            {
                if (di.codim() != cc)
                {
                    continue;
                }
                if (patchDomain.contains(patch + di))
                {
                    continue;
                }
                Point shrinkDirs;
                for (int d = 0; d < DIM; d++)
                {
                    shrinkDirs[d] = (di[d] != 0) ? -1 : 0;
                }
                Box k = K.grow(shrinkDirs).shift(di);
                PROTO_ASSERT(k.contains(di),
                             "MBLevelBoxData::define | Error: Data corruption");
                auto bounds = a_layout.domain().boundaries(block, di);
                Box boundaryPatches = patchDomain.adjacent(di, 1);
                for (auto bound : bounds)
                {
                    PROTO_ASSERT(bound.srcToDst == di,
                                 "MBLevelBoxData::define | Error: Data corruption");
                    auto adjBlock = bound.dstBlock;
                    Point adjDir = bound.dstToSrc;
                    auto adjLocalLayout = a_layout.blockLayout(adjBlock);
                    Box adjPatchDomain = adjLocalLayout.patchDomain().box();
                    Box adjBoundaryPatches = adjPatchDomain.edge(adjDir);
                    auto R = a_layout.domain().graph().rotation(block, di, adjBlock);
                    for (auto ki : k)
                    {
                        boundsKey_t key(patch, block, ki);
                        Point neighbor = patch + ki;
                        if (boundaryPatches.contains(neighbor))
                        {
                            auto adjPatch = R.rotateCell(neighbor, boundaryPatches, adjBoundaryPatches);
                            auto adjIndex = adjLocalLayout.find(adjPatch);

                            // for codim > 1, the rotation isn't defined in the directions
                            // that are linearly independent from the boundary surface
                            // to get this right, we rotate the in-plane coordinates and
                            // use the reverse arc for the out of plane coordinates

                            Point kiAdj = a_layout.domain().graph().reverseDir(block, adjBlock, ki);

                            // TODO: make this work with anisotropic ghost regions
                            Box patchBoundary = patchBox.adjacent(ki, ghostSize);
                            Box adjPatchBoundary = adjLocalLayout[adjIndex].edge(kiAdj, ghostSize);

                            auto adjIter = a_layout.find(adjPatch, adjBlock);
                            PROTO_ASSERT(adjIter != *a_layout.end(),
                                         "MBLevelBoxData::define | Error: Attempting to create \
                                    BoundaryData with nonexistent patch");
                            MBBoundaryData<T, C, MEM> boundData(
                                iter, adjIter, patchBoundary, adjPatchBoundary,
                                R.inverse(), a_boundGhost);
                            m_bounds[key].push_back(boundData);
                        }
                    }
                }
            }
        }
    }
    m_exchangeCopier.define(MBLevelExchangeCopierOp<T, C, MEM, CTR>(*this));
}

// HELPER FUNCTION FOR INITIALIZE
template<typename IN>
const IN&
parse_mb_level_arg(
        MBIndex& a_index, 
        Centering C,
        const IN& a_arg)
{
    return a_arg;
}

// HELPER FUNCTION FOR INITIALIZE
template<typename T, unsigned int C, MemType MEM, Centering CTR>
const BoxData<T,C,MEM>&
parse_mb_level_arg(
        MBIndex& a_index, 
        Centering _CTR,
        const MBLevelBoxData<T,C,MEM,CTR>& a_arg)
{
    return a_arg[a_index];
}

// HELPER FUNCTION FOR INITIALIZE
template<typename Func, MemType MEM>
const BoxData<double,DIM,MEM>
parse_mb_level_arg(
        MBIndex& a_index, 
        Centering C,
        const MBMap<Func, MEM>& a_arg)
{
    auto& layout = a_arg.map().layout();
    Box B0 = layout[a_index];
    B0 = B0.grow(a_arg.map().ghost()[0]);
    B0 = B0.grow(C);
    if (C == PR_NODE) {
        BoxData<double, DIM, MEM> data(B0);
        a_arg.map()[a_index].copyTo(data);
        return data;
    } else {
        auto block = layout.block(a_index);
        return a_arg(B0, block, C);
    }
}

// HELPER FUNCTION FOR INITIALIZE
template<int I = 0, typename Func, typename... LArgs, typename... FArgs>
typename std::enable_if<I == sizeof...(LArgs), void>::type
call_mb_level_forall(Func& a_func, MBIndex& a_index, Centering CTR,
        std::tuple<LArgs...> a_args,
        FArgs&&... a_fargs)
{
    forallInPlace_p(a_func, a_fargs...); 
}

// HELPER FUNCTION FOR INITIALIZE
template<int I = 0, typename Func, typename... LArgs, typename... FArgs>
typename std::enable_if<I < sizeof...(LArgs), void>::type
call_mb_level_forall(Func& a_func, MBIndex& a_index, Centering CTR,
        std::tuple<LArgs...> a_args,
        FArgs&&... a_fargs)
{
    auto& arg = parse_mb_level_arg(a_index, CTR, std::get<I>(a_args));
    call_mb_level_forall<I+1>(a_func, a_index, CTR, a_args, a_fargs..., arg);
}

template <typename T, unsigned int C, MemType MEM, Centering CTR>
template <typename Func, typename... Srcs>
void MBLevelBoxData<T, C, MEM, CTR>::initialize(Func &a_func, Srcs&... a_srcs)
{
    auto srcs = std::tuple<Srcs&...>(a_srcs...);
    for (auto iter : layout())
    {
        auto& patch = (*this)[iter];
        auto block = layout().block(iter);
        call_mb_level_forall(a_func, iter, CTR, srcs, patch, block);
    }
    /*
    for (int bi = 0; bi < numBlocks(); bi++)
    {
        auto &data = blockData(bi);
        data.initialize(a_func, bi, a_srcs...);
    }
    */
}

template <typename T, unsigned int C, MemType MEM, Centering CTR>
template <typename Func, typename... Srcs>
void MBLevelBoxData<T, C, MEM, CTR>::initConvolve(Func &a_func, Srcs... a_srcs)
{
    for (int bi = 0; bi < numBlocks(); bi++)
    {
        auto &data = blockData(bi);
        data.initConvolve(a_func, bi, a_srcs...);
    }
}

template <typename T, unsigned int C, MemType MEM, Centering CTR>
void MBLevelBoxData<T, C, MEM, CTR>::setVal(T a_value)
{
    for (auto data : m_data)
    {
        data->setVal(a_value);
    }
    fillBoundaries();
}

template <typename T, unsigned int C, MemType MEM, Centering CTR>
void MBLevelBoxData<T, C, MEM, CTR>::fillBoundaries()
{
    PR_TIME("MBLevelBoxData::fillBoundaries");
#ifdef PR_MPI
    {
        PR_TIME("MBLevelBoxData::fillBoundaries (MPI barrier)");
        barrier();
    }
#endif
    m_exchangeCopier.execute();
}

template <typename T, unsigned int C, MemType MEM, Centering CTR>
void MBLevelBoxData<T, C, MEM, CTR>::interpolate()
{
}

template <typename T, unsigned int C, MemType MEM, Centering CTR>
void MBLevelBoxData<T, C, MEM, CTR>::exchange()
{
    fillBoundaries();
    interpolate();
}

// TODO: This may be too much overhead in which case MBLevelBoxData
//   will likely require it's own custom Copier for copyTo
template <typename T, unsigned int C, MemType SRC_MEM, Centering CTR>
template <MemType DST_MEM>
void MBLevelBoxData<T, C, SRC_MEM, CTR>::copyTo(MBLevelBoxData<T, C, DST_MEM, CTR> &a_rhs) const
{
    PR_TIME("MBLevelBoxData::copyTo");
    PROTO_ASSERT(a_rhs.layout().domain().graph() == layout().domain().graph(),
                 "MBLevelBoxData::copyTo | Error: Incompatible datasets");
    /*
    for (int bi = 0; bi < numBlocks(); bi++)
    {
        (*this)[bi].copyTo(a_rhs[bi]);
    }
    */
    if (SRC_MEM == DST_MEM)
    {
        auto mutableSrc = const_cast<MBLevelBoxData<T, C, SRC_MEM, CTR> *>(this);
        MBLevelCopierOp<T, C, SRC_MEM, DST_MEM, CTR> op(*mutableSrc, a_rhs);
        MBLevelCopier<T, C, SRC_MEM, DST_MEM, CTR> copier;
        copier.define(op);
        copier.execute();
    }
    else
    {
        MayDay<void>::Abort("MBLevelBoxData::copyTo is not implemented for SRC_MEM != DST_MEM");
        // MBLevelBoxData<T, C, DST_MEM, CTR> tmpSrc(this->layout(), this->ghost());
        // this->copyToSimple(tmpSrc);
        // protoThreadSynchronize();
        // tmpSrc.copyTo(a_dest);
    }
}

template <typename T, unsigned int C, MemType MEM, Centering CTR>
BoxData<T, C, MEM> &
MBLevelBoxData<T, C, MEM, CTR>::operator[](MBIndex a_index)
{
    unsigned int block = m_layout->block(a_index);
    auto localIndex = m_layout->blockIndex(a_index);
    return (*m_data[block])[localIndex];
}

template <typename T, unsigned int C, MemType MEM, Centering CTR>
const BoxData<T, C, MEM> &
MBLevelBoxData<T, C, MEM, CTR>::operator[](MBIndex a_index) const
{
    unsigned int block = m_layout->block(a_index);
    auto localIndex = m_layout->blockIndex(a_index);
    return (*m_data[block])[localIndex];
}

template <typename T, unsigned int C, MemType MEM, Centering CTR>
Array<Point, DIM + 1>
MBLevelBoxData<T, C, MEM, CTR>::ghost() const
{
    return m_ghost;
}

template <typename T, unsigned int C, MemType MEM, Centering CTR>
LevelBoxData<T, C, MEM, CTR> &
MBLevelBoxData<T, C, MEM, CTR>::blockData(unsigned int a_block)
{
    return *m_data[a_block];
}

template <typename T, unsigned int C, MemType MEM, Centering CTR>
const LevelBoxData<T, C, MEM, CTR> &
MBLevelBoxData<T, C, MEM, CTR>::blockData(unsigned int a_block) const
{
    return *m_data[a_block];
}

template <typename T, unsigned int C, MemType MEM, Centering CTR>
unsigned int
MBLevelBoxData<T, C, MEM, CTR>::patchSize(unsigned int a_block) const
{
    return blockData(a_block).patchSize();
}

template <typename T, unsigned int C, MemType MEM, Centering CTR>
Box MBLevelBoxData<T, C, MEM, CTR>::patchBox(const MBIndex &a_index) const
{
    unsigned int block = layout().block(a_index);
    auto blockIndex = layout().blockIndex(a_index);
    return blockData(block).patchBox(blockIndex);
}

template <typename T, unsigned int C, MemType MEM, Centering CTR>
size_t
MBLevelBoxData<T, C, MEM, CTR>::linearSize() const
{
    size_t ret = 0;
    for (auto iter : layout())
    {
        ret += linearSize(iter);
    }
    return ret;
}

template <typename T, unsigned int C, MemType MEM, Centering CTR>
size_t
MBLevelBoxData<T, C, MEM, CTR>::linearSize(
    const MBIndex &a_index) const
{
    auto block = layout().block(a_index);
    return patchSize(block) * sizeof(T);
}

template <typename T, unsigned int C, MemType MEM, Centering CTR>
void MBLevelBoxData<T, C, MEM, CTR>::linearOut(
    void *a_buffer,
    const MBIndex &a_index,
    const Box &a_box) const
{
    PR_TIME("MBLevelBoxData::linearOut");
    const auto &patch = (*this)[a_index];
    CInterval comps(0, C - 1);
    patch.linearOut(a_buffer, a_box, comps);
}

template <typename T, unsigned int C, MemType MEM, Centering CTR>
void MBLevelBoxData<T, C, MEM, CTR>::linearOut(
    void *a_buffer,
    const MBIndex &a_index) const
{
    linearOut(a_buffer, a_index, (*this)[a_index].box());
}

template <typename T, unsigned int C, MemType MEM, Centering CTR>
void MBLevelBoxData<T, C, MEM, CTR>::linearOut(void *a_buffer) const
{
    char *writePoint = (char *)a_buffer;
    for (auto iter : layout())
    {
        linearOut(writePoint, iter);
        writePoint += linearSize(iter);
    }
}

template <typename T, unsigned int C, MemType MEM, Centering CTR>
size_t
MBLevelBoxData<T, C, MEM, CTR>::offset(unsigned int a_proc) const
{
    PR_TIME("MBLevelBoxData::offset");
    const auto &partition = layout().partition().partition();
    size_t ret = 0;
    unsigned int proc = 0;
    while (proc < a_proc)
    {
        for (auto pdata : partition)
        {
            if (pdata.second == proc)
            {
                auto block = pdata.first.second;
                ret += blockData(block).patchSize();
            }
        }
        proc++;
    }
    return ret;
}

template <typename T, unsigned int C, MemType MEM, Centering CTR>
std::vector<MBBoundaryData<T, C, MEM>> &
MBLevelBoxData<T, C, MEM, CTR>::bounds(MBIndex a_index, Point a_dir)
{
    auto patch = m_layout->point(a_index);
    auto block = m_layout->block(a_index);
    return m_bounds[boundsKey_t(patch, block, a_dir)];
}

template <typename T, unsigned int C, MemType MEM, Centering CTR>
const std::vector<MBBoundaryData<T, C, MEM>> &
MBLevelBoxData<T, C, MEM, CTR>::bounds(MBIndex a_index, Point a_dir) const
{
    auto patch = m_layout->point(a_index);
    auto block = m_layout->block(a_index);
    return m_bounds[boundsKey_t(patch, block, a_dir)];
}

template <typename T, unsigned int C, MemType MEM, Centering CTR>
bool MBLevelBoxData<T, C, MEM, CTR>::isBlockBoundary(
    MBIndex a_index, Point a_dir, unsigned int a_block)
{
    auto &boundaries = bounds(a_index, a_dir);
    for (auto b : boundaries)
    {
        if (m_layout->block(b.adjIndex) == a_block)
        {
            return true;
        }
    }
    return false;
}

template <typename T, unsigned int C, MemType MEM, Centering CTR>
const MBBoundaryData<T, C, MEM> &
MBLevelBoxData<T, C, MEM, CTR>::bound(MBIndex a_index, Point a_dir, unsigned int a_block) const
{
    PROTO_ASSERT(isBlockBoundary(a_index, a_dir, a_block),
                 "MBLevelBoxData::bound | Error: Inputs do not correspond to a block boundary.");
    auto &boundaries = bounds(a_index, a_dir);
    // for (auto b : boundaries)
    for (int bi = 0; bi < boundaries.size(); bi++)
    {
        auto &b = boundaries[bi];
        // END TEST CODE
        if (m_layout->block(b.adjIndex) == a_block)
        {
            return b;
        }
    }
    MayDay<void>::Error("MBLevelBoxData::bound | Error: Data corruption.");
}

template <typename T, unsigned int C, MemType MEM, Centering CTR>
MBBoundaryData<T, C, MEM> &
MBLevelBoxData<T, C, MEM, CTR>::bound(MBIndex a_index, Point a_dir, unsigned int a_block)
{
    PROTO_ASSERT(isBlockBoundary(a_index, a_dir, a_block),
                 "MBLevelBoxData::bound | Error: Inputs do not correspond to a block boundary.");
    auto &boundaries = bounds(a_index, a_dir);
    // for (auto b : boundaries)
    for (int bi = 0; bi < boundaries.size(); bi++)
    {
        auto &b = boundaries[bi];
        if (m_layout->block(b.adjIndex) == a_block)
        {
            return b;
        }
    }
    MayDay<void>::Error("MBLevelBoxData::bound | Error: Data corruption.");
}

template <typename T, unsigned int C, MemType MEM, Centering CTR>
bool
MBLevelBoxData<T, C, MEM, CTR>::compatible(const MBDisjointBoxLayout& a_layout) const
{
    return (m_layout->compatible(a_layout));
}

template <typename T, unsigned int C, MemType MEM, Centering CTR>
Var<T, C, MEM>
MBLevelBoxData<T, C, MEM, CTR>::operator[](const MBDataPoint &a_point)
{
    if (!a_point.inBoundary())
    {
        // bulk data is requested
        auto &patch = (*this)[a_point.index];
        PROTO_ASSERT(patch.box().contains(a_point.point),
                     "MBLevelBoxData::operator[](MBDataPoint) | Error: \
                The requested MBDataPoint is not valid.");
        return patch.var(a_point.point);
    }
    else
    {
        PROTO_ASSERT(a_point.block >= 0 && a_point.block < numBlocks(),
                     "MBLevelBoxData::operator[](MBDataPoint) | Error: \
                invalid boundary block specified: %i",
                     a_point.block);
        auto &boundData = (*this).bound(
            a_point.index, a_point.boundaryDir, a_point.block);
        PROTO_ASSERT(boundData.localData->box().contains(a_point.point),
                     "MBLevelBoxData::operator[](MBDataPoint) | Error: \
                The requested MBDataPoint is not valid.");
        return boundData.localData->var(a_point.point);
    }
}

template <typename T, unsigned int C, MemType MEM, Centering CTR>
Var<T, C, MEM>
MBLevelBoxData<T, C, MEM, CTR>::operator[](const MBDataPoint &a_point) const
{
    MBLevelBoxData<T, C, MEM, CTR> *ptr = const_cast<MBLevelBoxData<T, C, MEM, CTR> *>(*this);
    return (*ptr)[a_point];
}

template <typename T, unsigned int C, MemType MEM, Centering CTR>
std::vector<MBDataPoint>
MBLevelBoxData<T, C, MEM, CTR>::interpFootprint(
    const Point &a_center,
    const Point& a_ghost,
    const std::vector<Point> &a_footprint,
    const MBIndex &a_index)
{
    auto &srcPatch = (*this)[a_index];
    Box patchDomain = this->layout()[a_index].grow(a_ghost);

    PROTO_ASSERT(patchDomain.contains(a_center),
            "MBLevelBoxData::interpFootprint | Error: \
            Invalid center Point");
    auto block = layout().block(a_index);
    Box blockDomain = layout().blockLayout(block).domain().box();
    std::set<MBDataPoint> footprint;
    for (auto shift : a_footprint)
    {
        Point p = a_center + shift;
        if (blockDomain.contains(p) && patchDomain.contains(p))
        {
            MBDataPoint data(a_index, p);
            footprint.insert(data);
        }
        else if (!blockDomain.contains(p))
        {
            for (auto dir : Box::Kernel(1))
            {
                if (dir == Point::Zeros()) { continue; }
                auto boundData = bounds(a_index, dir);
                for (auto bi : boundData)
                {
                    if (bi.localBox().contains(p))
                    {
                        MBDataPoint data(a_index, p, dir, layout().block(bi.adjIndex));
                        footprint.insert(data);
                    }
                }
            }
        } else {
            // skip this point. It's in an adjacent patch that may not be on this proc
        }

    }
    std::vector<MBDataPoint> footprint_v(footprint.begin(), footprint.end());
    return footprint_v;
}

template <typename T, unsigned int C, MemType MEM, Centering CTR>
void MBLevelBoxData<T, C, MEM, CTR>::printBounds(const MBIndex &a_index, int a_prec)
{
    Box dirs = Box::Kernel(1);
    for (auto dir : dirs)
    {
        auto boundData = bounds(a_index, dir);
        for (auto bi : boundData)
        {
            pout() << "Dir: " << dir << " | Block: " << layout().block(bi.adjIndex) << std::endl;
            bi.localData->printData(a_prec);
        }
    }
}

template <typename T, unsigned int C, MemType MEM, Centering CTR>
MBIterator
MBLevelBoxData<T, C, MEM, CTR>::begin()
{
    return m_layout->begin();
}

template <typename T, unsigned int C, MemType MEM, Centering CTR>
MBIterator
MBLevelBoxData<T, C, MEM, CTR>::end()
{
    return m_layout->end();
}

//=======================================================================================
/// LEVEL EXCHANGE COPIER OP

template <typename T, unsigned int C, MemType MEM, Centering CTR>
MBLevelExchangeCopierOp<T, C, MEM, CTR>::MBLevelExchangeCopierOp(
    MBLevelBoxData<T, C, MEM, CTR> &a_data)
{
    m_data = &a_data;
}

template <typename T, unsigned int C, MemType MEM, Centering CTR>
uint64_t
MBLevelExchangeCopierOp<T, C, MEM, CTR>::linearSize(
    const Box &a_box,
    const MBIndex &a_index) const
{
    // TODO: Does this do the right thing?
    return C * a_box.size() * sizeof(T);
}

template <typename T, unsigned int C, MemType MEM, Centering CTR>
void MBLevelExchangeCopierOp<T, C, MEM, CTR>::linearOut(
    void *a_buffer, const MBMotionItem &a_info) const
{
    MemType mem = pointerMemType(a_buffer);
    PROTO_ASSERT(mem == MEM,
                 "MBLevelExchangeCopierOp::linearOut | Error: Memtype mismatch. Detected %i but should be %i",
                 mem, MEM);
    auto &patch = (*m_data)[a_info.fromIndex];
    CInterval cint(0, C - 1);
    patch.linearOut(a_buffer, a_info.fromRegion, cint);
}

template <typename T, unsigned int C, MemType MEM, Centering CTR>
void MBLevelExchangeCopierOp<T, C, MEM, CTR>::linearIn(
    void *a_buffer, const MBMotionItem &a_info) const
{
    MemType mem = pointerMemType(a_buffer);
    PROTO_ASSERT(mem == MEM,
                 "MBLevelCopierOp::linearIn | Error: Memtype mismatch. Detected %i but should be %i",
                 mem, MEM);
    const auto &srcIndex = a_info.fromIndex;
    const auto &dstIndex = a_info.toIndex;
    Box dirs = Box::Kernel(1);
    CInterval cint(0, C - 1);
    for (auto dir : dirs)
    {
        auto bounds = m_data->bounds(dstIndex, dir);
        for (auto bound : bounds)
        {
            PROTO_ASSERT(bound.localIndex == dstIndex,
                         "MBLevelExchangeCopierOp::linearIn | Error: Data corruption.");
            if (bound.adjIndex == srcIndex)
            {
                auto &boundData = *bound.localData;
                Box adjBox = bound.adjData->box();
                Box localBox = bound.localData->box();
                auto &R = bound.adjToLocal;
                R.rotateBuffer((T *)a_buffer, boundData.data(), adjBox, localBox);
                return;
            }
        }
    }
    MayDay<void>::Warning("MBLevelExchangeCopierOp::linearIn | Warning: \
            could not find a boundary to copy data into");
}

template <typename T, unsigned int C, MemType MEM, Centering CTR>
void MBLevelExchangeCopierOp<T, C, MEM, CTR>::localCopy(
    const MBMotionItem &a_info) const
{
    const auto &domain = a_info.fromRegion;
    const auto &srcIndex = a_info.fromIndex;
    const auto &range = a_info.toRegion;
    const auto &dstIndex = a_info.toIndex;

    auto &layout = m_data->layout();
    auto srcBlock = layout.block(srcIndex);
    auto dstBlock = layout.block(dstIndex);
    PROTO_ASSERT(srcBlock != dstBlock,
                 "MBLevelExchangeCopier::localCopy | Error: \
            source and destination blocks are the same.");
    bool linearizationTest = false;
    if (!linearizationTest)
    {
        PR_TIMERS("MBLevelCopierOp::localCopy");
        auto &src = (*m_data)[srcIndex];
        Box dirs = Box::Kernel(1);
        for (auto dir : dirs)
        {
            auto bounds = m_data->bounds(dstIndex, dir);
            for (auto bi : bounds)
            {
                PROTO_ASSERT(bi.localIndex == dstIndex,
                             "MBLevelExchangeCopierOp::localCopy | Error: Data corruption.");

                if (bi.adjIndex == srcIndex)
                {
                    auto &tmp = *bi.adjData;
                    auto &dst = *bi.localData;
                    auto &R = bi.adjToLocal;
                    PROTO_ASSERT(src.box().contains(tmp.box()),
                                 "MBLevelExchangeCopierOp::localCopy | Error: Data corruption.");
                    // TODO: Write more clever code that can do this in a single copy
                    src.copyTo(tmp);
                    tmp.copyTo(dst, R);
                    return;
                }
            }
        }
        MayDay<void>::Warning("MBLevelExchangeCopierOp::localCopy | Warning: \
                could not find a boundary to copy data into");
    }
    else
    {
        PR_TIMERS("MBLevelExchangeCopierOp::localCopy (linearization test)");
        int sizeSrc = linearSize(domain, srcIndex);
        int sizeDst = linearSize(range, dstIndex);

        if (sizeSrc != sizeDst)
        {
            MayDay<void>::Abort(
                "LinearizationTest failure: dest and source have different sizes");
        }

        void *srcBuffer = proto_malloc<MEM>(sizeSrc);
        void *dstBuffer = proto_malloc<MEM>(sizeDst);
        {
            PR_TIMERS("linearIn / linearOut inside LDOP::op");
            linearOut(srcBuffer, a_info);
            proto_memcpy<MEM, MEM>(srcBuffer, dstBuffer, sizeSrc);
            linearIn(dstBuffer, a_info);
        }
    }
}

template <typename T, unsigned int C, MemType MEM, Centering CTR>
void MBLevelExchangeCopier<T, C, MEM, CTR>::buildMotionPlans(
    MBLevelExchangeCopierOp<T, C, MEM, CTR> &a_op)
{
    auto layout = a_op.m_data->layout();
    auto ghostArray = a_op.m_data->ghost();
    auto &graph = layout.domain().graph();
    Box dirs = Box::Kernel(1);
    int localRank = procID();
    for (auto iter : layout)
    {
        auto localBlock = layout.block(iter);
        for (auto dir : dirs)
        {
            auto bounds = a_op.m_data->bounds(iter, dir);
            for (auto bound : bounds)
            {
                PROTO_ASSERT(layout.procID(bound.localIndex) == localRank,
                             "MBLevelExchangeCopier::buildMotionPlans | Error: Data corruption.");
                int adjRank = layout.procID(bound.adjIndex);
                auto adjBlock = layout.block(bound.adjIndex);
                auto boundCodim = graph.connectivity(localBlock, adjBlock).codim();
                Point ghost = ghostArray[boundCodim];
                Point adjDir = graph.reverseDir(localBlock, adjBlock, dir);

                Box adjBox_to = bound.adjData->box();
                Box localBox_to = bound.localData->box();
                MBMotionItem toMotionItem(bound.adjIndex, bound.localIndex, adjBox_to, localBox_to);
                if (localRank == adjRank)
                {
                    this->m_localMotionPlan.push_back(toMotionItem);
                }
                else
                {
                    toMotionItem.procID = adjRank;
                    this->m_toMotionPlan.push_back(toMotionItem);

                    Box adjBox_from = adjBox_to.shift(adjDir * ghost);
                    Box localBox_from = localBox_to.shift(-dir * ghost);
                    MBMotionItem fromMotionItem(bound.localIndex, bound.adjIndex, localBox_from, adjBox_from);
                    fromMotionItem.procID = adjRank;
                    this->m_fromMotionPlan.push_back(fromMotionItem);
                }
            }
        }
    }
    this->sort();
}

template <typename T, unsigned int C, MemType SRC_MEM, MemType DST_MEM, Centering CTR>
MBLevelCopierOp<T, C, SRC_MEM, DST_MEM, CTR>::MBLevelCopierOp(
    MBLevelBoxData<T, C, SRC_MEM, CTR> &a_src,
    MBLevelBoxData<T, C, DST_MEM, CTR> &a_dst)
{
    PROTO_ASSERT(SRC_MEM == DST_MEM,
                 "MBLevelCopierOp::constructor | Error: \
        MBLevelCopierOp with different MemType src and dst is not yet implemented.");
    PROTO_ASSERT(a_src.layout().domain() == a_dst.layout().domain(),
                 "MBLevelCopier::Constructor | Error: \
            src and dst must have the same problem domain");
    m_src = &a_src;
    m_dst = &a_dst;
}

template <typename T, unsigned int C, MemType SRC_MEM, MemType DST_MEM, Centering CTR>
uint64_t
MBLevelCopierOp<T, C, SRC_MEM, DST_MEM, CTR>::linearSize(
    const Box &a_box,
    const MBIndex a_index) const
{
    // TODO: Does this do the right thing?
    return C * a_box.size() * sizeof(T);
}

template <typename T, unsigned int C, MemType SRC_MEM, MemType DST_MEM, Centering CTR>
void MBLevelCopierOp<T, C, SRC_MEM, DST_MEM, CTR>::linearOut(
    void *a_buffer, const MBMotionItem &a_info) const
{
    MemType mem = pointerMemType(a_buffer);
    PROTO_ASSERT(mem == SRC_MEM,
                 "MBLevelCopierOp::linearOut | Error: Memtype mismatch. Detected %i but should be %i",
                 mem, SRC_MEM);
    auto &patch = (*m_src)[a_info.fromIndex];
    CInterval cint(0, C - 1);
    patch.linearOut(a_buffer, a_info.fromRegion, cint);
}

template <typename T, unsigned int C, MemType SRC_MEM, MemType DST_MEM, Centering CTR>
void MBLevelCopierOp<T, C, SRC_MEM, DST_MEM, CTR>::linearIn(
    void *a_buffer, const MBMotionItem &a_info) const
{
    MemType mem = pointerMemType(a_buffer);
    PROTO_ASSERT(mem == DST_MEM,
                 "MBLevelCopierOp::linearIn | Error: Memtype mismatch. Detected %i but should be %i",
                 mem, DST_MEM);
    auto &patch = (*m_dst)[a_info.toIndex];
    CInterval cint(0, C - 1);
    patch.linearIn(a_buffer, a_info.toRegion, cint);
}

template <typename T, unsigned int C, MemType SRC_MEM, MemType DST_MEM, Centering CTR>
void MBLevelCopierOp<T, C, SRC_MEM, DST_MEM, CTR>::localCopy(
    const MBMotionItem &a_info) const
{
    bool linearizationTest = false;

    const auto &domain = a_info.fromRegion;
    const auto &srcIndex = a_info.fromIndex;
    const auto &range = a_info.toRegion;
    const auto &dstIndex = a_info.toIndex;

    if (!linearizationTest)
    {
        PR_TIMERS("MBLevelCopierOp::localCopy");
        auto &src = (*m_src)[srcIndex];
        auto &dst = (*m_dst)[dstIndex];
        Point shift = range.low() - domain.low();
        src.copyTo(dst, domain, shift);
    }
    else
    {
        PR_TIMERS("MBLevelCopierOp::localCopy (linearization test)");
        int sizeSrc = linearSize(domain, srcIndex);
        int sizeDst = linearSize(range, dstIndex);

        if (sizeSrc != sizeDst)
        {
            MayDay<void>::Abort(
                "LinearizationTest failure: dest and source have different sizes");
        }

        void *srcBuffer = proto_malloc<SRC_MEM>(sizeSrc);
        void *dstBuffer = proto_malloc<DST_MEM>(sizeDst);
        {
            PR_TIMERS("linearIn / linearOut inside LDOP::op");
            linearOut(srcBuffer, a_info);
            proto_memcpy<SRC_MEM, DST_MEM>(srcBuffer, dstBuffer, sizeSrc);
            if (SRC_MEM != DST_MEM)
            {
#ifdef PROTO_ACCEL
                protoThreadSynchronize();
#endif
            }
            linearIn(dstBuffer, a_info);
        }
    }
}

template <typename T, unsigned int C, MemType SRC_MEM, MemType DST_MEM, Centering CTR>
void MBLevelCopier<T, C, SRC_MEM, DST_MEM, CTR>::buildMotionPlans(
    MBLevelCopierOp<T, C, SRC_MEM, DST_MEM, CTR> &a_op)
{
    PROTO_ASSERT(SRC_MEM == DST_MEM,
                 "MBLevelCopier::buildMotionPlans | Error: \
        Copy operations between HOST and DEVICE are not yet implemented.");
    auto srcLayout = a_op.m_src->layout();
    auto dstLayout = a_op.m_dst->layout();
    auto dstGhost = a_op.m_dst->ghost()[0];
    // For each Box on this processor, compute all the pairs of motionItems that move data
    // to this Box.
    for (auto diter : dstLayout)
    {
        unsigned int block = dstLayout.block(diter);
        auto srcBlockLayout = srcLayout.blockLayout(block);
        auto dstBlockLayout = dstLayout.blockLayout(block);
        auto srcBlockDomain = srcBlockLayout.domain();
        auto dstBlockDomain = dstBlockLayout.domain();
        auto srcBoxSize = srcBlockLayout.boxSize();
        auto dstBoxSize = dstBlockLayout.boxSize();

        // destGhost[0] is the ghost size in the bulk of each block
        // TODO: Make this anisotropic
        Box bx = dstLayout[diter].grow(dstGhost);

        // Compute all the patches (as Points in a bitmap) that might move data to bx.

        Box srcPatches = bx.coarsen(srcBoxSize);
        // Iterate over the possible source patches to contruct MotionItems.

        for (auto pt : srcPatches)
        {
            // Find the possible source patch, or its periodic image.

            auto di = srcLayout.find(pt, block);
            if (di != *srcLayout.end())
            {
                Box destBx = Box(pt, pt).refine(srcBoxSize) & bx;
                if (!destBx.empty())
                {
                    // Create MotionItem.
                    // Compute srcBx to be the shifted image of destBx. If destBx
                    // is in the ProblemDomain box, srcBx == destBx.
                    Box srcBx = Box(srcBlockDomain.image(destBx.low()),
                                    srcBlockDomain.image(destBx.high()));
                    MBMotionItem toMotionItem(di, diter, srcBx, destBx);
                    int srcProcID = srcLayout.procID(di);
                    if (srcProcID == procID())
                    {
                        PROTO_ASSERT(srcProcID == dstLayout.procID(diter),
                                     "LevelCopier::buildMotionPlan | Error: Data corruption.");
                        this->m_localMotionPlan.push_back(toMotionItem);
                    }
                    else
                    {
                        PROTO_ASSERT(srcProcID != dstLayout.procID(diter),
                                     "LevelCopier::buildMotionPlan | Error: Data corruption.");
                        toMotionItem.procID = srcProcID;
                        this->m_toMotionPlan.push_back(toMotionItem);
                    }
                }
            }
        }
    }

    // Iterate over Boxes in the source to construct fromMotionPlan
    for (auto siter : srcLayout)
    {
        unsigned int block = srcLayout.block(siter);
        auto srcBlockLayout = srcLayout.blockLayout(block);
        auto dstBlockLayout = dstLayout.blockLayout(block);
        auto srcBlockDomain = srcBlockLayout.domain();
        auto dstBlockDomain = dstBlockLayout.domain();
        auto srcBoxSize = srcBlockLayout.boxSize();
        auto dstBoxSize = dstBlockLayout.boxSize();

        // We grow the source Box by dstGhost to capture all the patches that the
        // source will copy to. dProbC contains all the points in the bitmap for
        // the destination.
        Box srcBox = srcLayout[siter];
        Box dProbC = srcBox.grow(dstGhost).coarsen(dstBoxSize);
        for (auto pt : dProbC)
        {
            auto di = dstLayout.find(pt, block);
            if (di != *dstLayout.end())
            {
                Box bxDestShifted = Box(pt, pt).refine(dstBoxSize).grow(dstGhost);
                Box bxDest = dstLayout[di].grow(dstGhost);
                Point shiftSrcDest = bxDest.low() - bxDestShifted.low();
                bxDest = bxDest & (srcBox.shift(shiftSrcDest));
                Box bxSrc = bxDestShifted & srcBox;
                if (!bxSrc.empty())
                {
                    int dstProcID = dstLayout.procID(di);
                    if (dstProcID != procID())
                    {
                        MBMotionItem fromMotionItem(siter, di, bxSrc, bxDest);
                        fromMotionItem.procID = dstProcID;
                        this->m_fromMotionPlan.push_back(fromMotionItem);
                    }
                }
            }
        }
    } // end src loop
    this->sort();
}

MBDataPoint::MBDataPoint(MBIndex a_index, Point a_point)
{
    index = a_index;
    point = a_point;
    block = -1;
    boundaryDir = Point::Zeros();
    m_inBoundary = false;
}
MBDataPoint::MBDataPoint(MBIndex a_index, Point a_point, Point a_boundaryDir, int a_block)
{
    index = a_index;
    block = a_block;
    boundaryDir = a_boundaryDir;
    point = a_point;
    m_inBoundary = true;
}

bool MBDataPoint::operator<(const MBDataPoint &a_rhs) const
{
    if (index.global() != a_rhs.index.global())
    {
        return index.global() < a_rhs.index.global();
    }
    else if (point != a_rhs.point)
    {
        return point < a_rhs.point;
    }
    if (m_inBoundary)
    {
        if (block != a_rhs.block)
        {
            return block < a_rhs.block;
        }
        else if (boundaryDir != a_rhs.boundaryDir)
        {
            return boundaryDir < a_rhs.boundaryDir;
        } else {
            return false;
        }
    }
    return false;
}

bool MBDataPoint::operator==(const MBDataPoint &a_rhs) const
{
    if (index != a_rhs.index)
    {
        return false;
    }
    if (point != a_rhs.point)
    {
        return false;
    }
    if (m_inBoundary)
    {
        if (block != a_rhs.block)
        {
            return false;
        }
        if (boundaryDir != a_rhs.boundaryDir)
        {
            return false;
        }
    }
    return true;
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
BoxData<T, C, MEM>& MBDataPoint::patch(MBLevelBoxData<T, C, MEM, CTR>& a_data) const
{
    if (m_inBoundary)
    {
        auto& boundData = a_data.bound(index, boundaryDir, block);
        return *boundData.localData;
    } else {
        return a_data[index];
    }
}
