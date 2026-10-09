MBInterpOp::MBInterpOp()
{
    PR_TIME("MBInterpOp::constructor");
}

MBInterpOp::MBInterpOp(Point a_ghost, unsigned int a_order)
{
    PR_TIME("MBInterpOp::constructor");
    m_order = a_order;
    m_ghost = a_ghost;
}

template<typename MAP, MemType MEM>
MBInterpOp::MBInterpOp(
        const MBLevelMap<MAP,MEM>&    a_map,
        int a_order)
{
    PR_TIME("MBInterpOp::constructor");
    define(a_map, a_order);
}

template<typename MAP, MemType MEM>
void MBInterpOp::define(
        const MBLevelMap<MAP,MEM>&    a_map,
        int a_order)
{
    MBInterpLayout interpLayout(a_map.layout(), a_order);
    define(a_map, interpLayout, a_order);
}

template<typename MAP, MemType MEM>
void MBInterpOp::define(
        const MBLevelMap<MAP,MEM>&  a_map,
        MBInterpLayout              a_interpLayout,
        int                         a_order)
{
    PR_TIME("MBInterpOp::define");
    PROTO_ASSERT(a_interpLayout.layout().compatible(a_map.layout()),
        "MBInterpOp::define | Error: Incompatible layouts among inputs");

    m_ops.clear();
    m_radialRuns.clear();
    m_points.clear();
    m_cubedSphere = false;
    m_order = a_order;
    m_ghost = a_map.ghost()[0];
    m_interpLayout = a_interpLayout;
    for (auto iter : m_interpLayout.layout())
    {
        for (auto dir : Box::Kernel(1))
        {
           addPointOps(iter, dir, a_map);
        }
    }
}

template<typename MAP, MemType MEM>
void MBInterpOp::define(
        const std::vector<std::shared_ptr<MBLevelMap<MAP,MEM>>>&  a_maps,
        MBInterpLayout    a_interpLayout,
        int                         a_order)
{
    PR_TIME("MBInterpOp::define");
    PROTO_ASSERT(a_maps.size() == a_interpLayout.layout().numBlocks(),
        "MBInterpOp::define | Error: Incorrect number of maps. There must be one map per block");
    for (auto bi = 0; bi < a_interpLayout.layout().numBlocks(); bi++)
    {
        PROTO_ASSERT(a_interpLayout.layout().compatible(a_maps[bi]->layout()),
            "MBInterpOp::define | Error: Incompatible layouts among inputs");
        PROTO_ASSERT(a_maps[bi]->ghost() == a_maps[0]->ghost(),
            "MBInterpOp::define | Error: Input maps have different ghost regions");
    }

    m_ops.clear();
    m_radialRuns.clear();
    m_points.clear();
    m_cubedSphere = false;
    m_order = a_order;
    m_ghost = a_maps[0]->ghost()[0];
    m_interpLayout = a_interpLayout;
    for (auto iter : m_interpLayout.layout())
    {
        auto block = m_interpLayout.layout().block(iter);
        auto map = a_maps[block];
        for (auto dir : Box::Kernel(1))
        {
           addPointOps(iter, dir, *map);
        }
    }
}

template<typename MAP, MemType MEM>
void MBInterpOp::addPointOps(const MBIndex& index, Point dir, const MBLevelMap<MAP,MEM>& map)
{
    PR_TIME("MBInterpOp::addPointOps");
    auto& layout = m_interpLayout.layout();
    auto patch = layout.patch(index);
    PatchID adjPatch = patch.point + dir;
    BlockIndex block = patch.block;
    auto& patchDomain = layout.patchDomain();

    if (patchDomain.isPointInInterior(adjPatch, block)) { return; }
    else if (patchDomain.isPointInDomainBoundary(adjPatch, block)) { return; }
    else if (layout.isRefinementBoundary(patch, dir)) { return; }
    else {
        Box patchBox = layout[index];
        Point targetGhost = m_ghost;
        if (m_cubedSphere) { targetGhost -= Point::Ones(2) - Point::Basis(0, 2); }
        Box boundBox = patchBox.adjacent(targetGhost * dir);
        for (auto boundPoint : boundBox)
        {
            if (this->containsOp(boundPoint, index))
            {
                continue;
            }
            //auto footprint = m_interpLayout.footprint(boundPoint, dir, index);
            addPointOp(boundPoint, index, map);
        }
    }
}

bool MBInterpOp::containsOp(Point dataPoint, const MBIndex& index) const
{
    PR_TIME("MBInterpOp::containsOp");
    MBPointID key(dataPoint, index.global());
    return (m_points.count(key) > 0);
}

template<typename MAP, MemType MEM>
void MBInterpOp::addPointOp(
    Point dataPoint,
    const MBIndex& index,
    const MBLevelMap<MAP,MEM>& map)
{
    PR_TIME("MBInterpOp::addPointOp");
    PROTO_ASSERT(!containsOp(dataPoint, index),
        "MBInterpOp::adPointOp | Error: Trying to add duplicate op");
    
    MBDataPoint dstDataPoint(index, dataPoint, m_interpLayout.layout());
    std::shared_ptr<MBPointInterpOp> op;
    if (m_cubedSphere)
    {
        if (map.layout().domain().getBlock(dstDataPoint.dstBlock()).box().containsPoint(dataPoint)) { return; }
        op = std::make_shared<MBPointInterpOp>(dstDataPoint, m_ghost, map, m_order);
    }
    else
    {
        op = std::make_shared<MBPointInterpOp>(dstDataPoint, m_ghost, map, m_interpLayout, m_order);
    }
    
    m_ops.push_back(op);
    MBPointID key(dataPoint, index.global());
    m_points.insert(key);

    copyAlongSymmetryAxes(*op);
}

void MBInterpOp::copyAlongSymmetryAxes(const MBPointInterpOp& pointOp)
{
    PR_TIME("MBInterpOp::copyAlongSymmetryAxis");
    Point dataPoint = pointOp.target().dstPoint();
    MBIndex index = pointOp.target().dstIndex();
    Box patchBox = m_interpLayout.layout()[index];
    Point boundaryDir = patchBox.whichBoundaryContains(dataPoint);
    Box patchBoundBox = patchBox.adjacent(boundaryDir * m_ghost);
    for (int dd = 0; dd < DIM; dd++)
    {
        PR_TIME("MBInterpOp::copyAlongSymmetryAxis");

        if (!m_interpLayout.isCopied(dd)) { continue; }

        Point low = dataPoint;
        Point high = dataPoint;
        low[dd] = patchBoundBox.low()[dd];
        high[dd] = patchBoundBox.high()[dd];
        Box axisBox(low, high);
        PROTO_ASSERT(axisBox.containsPoint(dataPoint),
            "MBInterpOp::copyAlongSymmetryAxes | Error: Data Corruption");
        for (auto axisPoint : axisBox)
        {
            if (containsOp(axisPoint, index)) { continue; }
            Point shift = axisPoint - dataPoint;
            addShiftedOp(pointOp, shift);
        }
    }
}

void MBInterpOp::addShiftedOp(
    const MBPointInterpOp& baseOp,
    Point shift)
{
    PR_TIME("MBInterpOp::addShiftedOp");
    PROTO_ASSERT(containsOp(baseOp.target().dstPoint(), baseOp.target().dstIndex()),
        "MBInterpOp::adPointOp | Error: Trying to add duplicate op");
    if (shift == Point::Zeros()) { return; }
    auto& layout = m_interpLayout.layout();
    MBDataPoint shiftDst = baseOp.target().shift(shift, layout);
    std::vector<MBDataPoint> shiftSrc;
    for (auto si : baseOp.sources())
    {
        shiftSrc.push_back(si.shift(shift, layout));
    }
    auto coefs = baseOp.coefs();
    auto shiftedOp = std::make_shared<MBPointInterpOp>(shiftDst, shiftSrc, coefs);
    
    m_ops.push_back(shiftedOp);
    MBPointID key(shiftDst.dstPoint(), baseOp.target().dstIndex().global());
    m_points.insert(key);
}

void MBInterpOp::insert(std::shared_ptr<MBPointInterpOp> a_op)
{
    PR_TIME("MBInterpOp::insert");
    m_radialRuns.clear();
    Point dst = a_op->target().dstPoint();
    int globalIndex = a_op->target().dstIndex();
    MBPointID key(dst, globalIndex);
    PROTO_ASSERT(m_points.count(key) == 0,
            "MBInterpOp::insert | Tried to insert a duplicate MBPointInterpOp");

    m_ops.push_back(a_op);
    m_points.insert(key);
}

template<typename T, unsigned int C, MemType MEM>
void MBInterpOp::apply(
        MBLevelBoxData<T, C, MEM>& a_dst,
        MBLevelBoxData<T, C, MEM>& a_src) const
{
    PR_TIME("MBInterpOp::apply");
    // Manually inserted operators can connect different boundary-data layouts.
    if (m_interpLayout.layout().defined())
    {
        PROTO_ASSERT(m_interpLayout.layout().compatible(a_src.layout()), "Incompatible Input Data");
        PROTO_ASSERT(m_interpLayout.layout().compatible(a_dst.layout()), "Incompatible Output Data");
    }
    if constexpr (MEM == HOST)
    {
        if (m_cubedSphere && !m_radialRuns.empty() && &a_dst != &a_src)
        {
            for (const auto& run : m_radialRuns)
            {
                const auto& op = m_ops[run.first];
                const Point first = op->target().dstPoint();
                const Box& box = a_dst[op->target().dstIndex()].box();
                Point probe = first;
                probe[0] = box.low()[0];
                if (!box.containsPoint(probe)) { continue; }
                const int low = std::max(first[0], box.low()[0]);
                const int high = std::min(first[0] + run.count - 1, box.high()[0]);
                if (low > high) { continue; }
                op->applyRadialRun(a_dst, a_src, low - first[0], high - low + 1);
            }
            return;
        }
    }
    for (auto& op : m_ops)
    {
        auto block = a_dst.layout().block(op->target().dstIndex());
        if (a_dst[op->target().dstIndex()].box().containsPoint(op->target().dstPoint()))
        {
            op->apply(a_dst, a_src);
        }
    }
}
void MBInterpOp::print() const
{
    pr_out() << "Number of PointInterpOps: " << m_ops.size() << std::endl;
    for (auto& op : m_ops)
    {
        op->print(); 
    }
}
void MBInterpOp::writeFootprint(std::string filename) const
{
    std::map<int, int> counters;
    for (auto& op : m_ops)
    {
        int block = op->target().dstBlock();
        int ii = counters[block]++;
        std::string suffix = ("B" + std::to_string(block) + "_I" + std::to_string(ii));
        op->writeFootprint(filename+suffix);
    }
}
template<typename MAP, MemType MEM>
void MBInterpOp::writeLevelFootprint(const MBLevelMap<MAP, MEM>& map, std::string filename) const
{
    std::map<int, int> counters;
    for (auto& op : m_ops)
    {
        int block = op->target().dstBlock();
        int ii = counters[block]++;
        std::string suffix = ("B" + std::to_string(block) + "_I" + std::to_string(ii));
        op->writeLevelFootprint(map, filename+suffix);
    }
}

template<typename T, unsigned int C, MemType MEM>
void MBInterpOp::printErrorPoints(
        MBLevelBoxData<T, C, MEM>& a_err,
        T                          a_eps)
{
    auto& layout = a_err.layout();
    for (auto iter : layout)
    {
        auto block = layout.block(iter);
        auto& err_i = a_err[iter];
        for (auto pi : err_i.box())
        {
            for (int cc = 0; cc < C; cc++)
            {
                if (err_i(pi, cc) > a_eps)
                {   
                    for (auto& op : m_ops)
                    {
                        auto& dst = op->target();
                        if (dst.dstPoint() == pi && dst.dstBlock() == block)
                        {
                            std::cout << "\nError Point found | error: " << err_i(pi, cc) << " | target: " << pi << " | block: " << block << std::endl;
                            std::cout << "Source Points: " << std::endl;
                            for (auto& si : op->sources())
                            {
                                std::cout << "\t" << si << std::endl;
                            }
                            break;
                        }
                    }
                    break;
                }
                
            }
        }
    }
}


template<typename MAP, typename T, unsigned int C, MemType MEM, Centering CTR>
void interpBoundaries(MBLevelBoxData<T, C, MEM, CTR>& a_data, unsigned int a_order)
{
    MBLevelMap<MAP, MEM> map(a_data.layout(), a_data.ghost());
    interpBoundaries<MAP, T, C, MEM, CTR>(a_data, map, a_order);
}

template<typename MAP, typename T, unsigned int C, MemType MEM, Centering CTR>
void interpBoundaries(
        MBLevelBoxData<T, C, MEM, CTR>& a_data,
        MBLevelMap<MAP, MEM>& a_map,
        unsigned int a_order)
{
    std::vector<Point> footprint;
    Box K = Box::Kernel(2);
    for (auto pi : K)
    {
        if (pi.abs().sum() <= 2)
        {
            footprint.push_back(pi);
        }
    }
    MBInterpOp op(a_data.ghost()[0], a_order);
    op.define(a_map, footprint);
    op.apply(a_data, a_data);
}







// Specialized construction is explicit so generic AMR interpolation keeps its
// own footprint selection and refinement-boundary handling.
template<typename MAP, MemType MEM>
void MBInterpOp::defineCubedSphere(
    const std::vector<std::shared_ptr<MBLevelMap<MAP, MEM>>>& maps,
    MBInterpLayout interpLayout, Point ghost, int order)
{
    PROTO_ASSERT(DIM == 3 && order == 4, "Cubed-sphere interpolation requires DIM=3 and order=4");
    PROTO_ASSERT(maps.size() == interpLayout.layout().numBlocks(), "Expected one map per block");
    m_ops.clear();
    m_radialRuns.clear();
    m_points.clear();
    m_cubedSphere = true;
    m_order = order;
    m_ghost = ghost;
    m_interpLayout = interpLayout;
    for (const auto& map : maps)
    {
        PROTO_ASSERT(map && map->layout().compatible(interpLayout.layout()), "Incompatible cubed-sphere map");
    }
    for (auto index : interpLayout.layout())
    {
        for (auto dir : Box::Kernel(1))
        {
            addPointOps(index, dir, *maps[interpLayout.layout().block(index)]);
        }
    }
    buildRadialRuns();
}

void MBInterpOp::buildRadialRuns()
{
    PR_TIME("MBInterpOp::buildRadialRuns");
    m_radialRuns.clear();
    const Point step = Point::Basis(0);
    auto shifted = [&](const MBDataPoint& a, const MBDataPoint& b)
    {
        return a.dstIndex() == b.dstIndex() && a.srcIndex() == b.srcIndex()
            && a.dstBlock() == b.dstBlock() && a.srcBlock() == b.srcBlock()
            && a.boundaryDir() == b.boundaryDir()
            && b.dstPoint() == a.dstPoint() + step
            && b.srcPoint() == a.srcPoint() + step;
    };
    for (size_t ii = 0; ii < m_ops.size(); ++ii)
    {
        bool append = ii > 0 && shifted(m_ops[ii-1]->target(), m_ops[ii]->target())
            && m_ops[ii-1]->coefs() == m_ops[ii]->coefs()
            && m_ops[ii-1]->sources().size() == m_ops[ii]->sources().size();
        if (append)
        {
            const auto& previous = m_ops[ii-1]->sources();
            const auto& current = m_ops[ii]->sources();
            for (size_t jj = 0; jj < current.size(); ++jj)
            {
                if (!shifted(previous[jj], current[jj])) { append = false; break; }
            }
        }
        if (append) { ++m_radialRuns.back().count; }
        else { m_radialRuns.push_back({ii, 1}); }
    }
}
