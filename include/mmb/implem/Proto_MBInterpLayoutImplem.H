
MBInterpLayout::MBInterpLayout(const MBDisjointBoxLayout& layout, int order)
{
    m_baseFootprint.clear();
    m_copyAxes.fill(false);
    m_layout = layout;
    m_baseFootprint.clear();
    switch (order)
    {
    case 4:
        for (auto pi : Box::Kernel(2))
        {
            if (pi.abs().sum() <= 2)
            {
                addPoint(pi);
            }
        }
        break;
    case 5:
        for (auto pi : Box::Kernel(2))
        {
            if (pi.abs().sum() <= 3)
            {
                addPoint(pi);
            }
        }
        break;
    default:
        MayDay<void>::Error("MBInterpLayout::define | Function undefined for order %i", order);
        break;
    }
}

void MBInterpLayout::ignoreAxis(int dir)
{
    std::vector<Point> footprint(m_baseFootprint.begin(), m_baseFootprint.end());
    m_baseFootprint.clear();
    for (auto fi : footprint)
    {
        if (fi[dir] == 0) { addPoint(fi); }
    }
}

std::set<Point> MBInterpLayout::footprint(
    Point center,
    MBIndex index) const
{
    auto patchDomain = m_layout.patchDomain();
    BlockIndex block = m_layout.block(index);
    PatchID patch = m_layout.point(index);
    Box patchBox = m_layout[index];
    Point boundaryDir = patchBox.whichBoundaryContains(center);
    PROTO_ASSERT(boundaryDir != Point::Zeros(), "Center point should be inside of the boundary of a patch");
    PatchID adjPatch = patch + boundaryDir;
    std::set<Point> footprint(m_baseFootprint.begin(), m_baseFootprint.end());
    addBoundaryFootprint(footprint, center, index);
    if (patchDomain.isPointInTriplePoint(adjPatch, block))
    {
        addTriplePointFootprint(footprint, center, boundaryDir, index);
    } else if (patchDomain.isPointInBlockBoundary(adjPatch, block))
    {
        for (auto orthDir : Point::DirectionsOfCodim(1))
        {
            bool isOrthDir = (orthDir.dot(boundaryDir) == 0);
            bool isTriplePointAdjacent = m_layout.domain().isPointInTriplePoint(center + orthDir, block);
            if (isTriplePointAdjacent && isOrthDir)
            {
                addTriplePointFootprint(footprint, center, boundaryDir + orthDir, index);
            }

        }
    }
    removeDomainBoundaryPoints(footprint, center, index);
    removeTriplePointPoints(footprint, center, index);
    return footprint;
}

template<typename MAP, MemType MEM>
inline std::vector<MBDataPoint>
MBInterpLayout::computeSourcePoints(
    std::set<Point> footprint,
    const MBLevelMap<MAP,MEM>& map,
    Point center,
    MBIndex index) const
{
    PR_TIME("MBInterpLayoutOp::computeSourcePoints");
    PROTO_ASSERT(map.layout().compatible(layout()),
        "MBInterpLayout::computeSourcePoints | Error: map layout is incompatible with InterpLayout internals");

    auto patch = layout().patch(index);
    Box patchBox = layout()[index];
    patchBox = patchBox.grow(map.ghost()[0]); // all valid data that will be locally available on the destination patch
    auto domain = layout().domain();
    BlockIndex block = patch.block;
    Box domainBox = layout().domain().getBlock(block).box();
    
    std::vector<MBDataPoint> srcs;
    std::set<Point> augmentedFootprint;
    for (auto fi : footprint)
    {
        Point pi = center + fi;
        bool success = addSourceFromShiftedPoint(pi, srcs, map, index);
        if (!success)
        {
            augmentFootprint(fi, augmentedFootprint, footprint);
        }
    }
    for (auto fi : augmentedFootprint)
    {
        Point pi = center + fi;
        addSourceFromShiftedPoint(pi, srcs, map, index);
    }
    PROTO_ASSERT(srcs.size() > 0,
        "MBInterpLayout::computeSourcePoints | Error: Could not compute any source data points for interpolation to point %s in block %i",
        center.str(), patch.block);
    return srcs;
}

void MBInterpLayout::augmentFootprint(
    Point shift,
    std::set<Point>& augmentedFootprint,
    const std::set<Point>& footprint) const
{
    Point f0 = -shift.dir();
    for (auto fi : footprint)
    {
        Point fj = fi + f0;
        if (footprint.count(fj) <= 0) { augmentedFootprint.insert(fj); }
    }
}

// IMPLEMENT ME
template<typename MAP, MemType MEM>
bool MBInterpLayout::addSourceFromShiftedPoint(
    Point shiftedPoint,
    std::vector<MBDataPoint>& sources,
    const MBLevelMap<MAP,MEM>& map,
    MBIndex index) const
{
    auto patch = layout().patch(index);
    Point pi = shiftedPoint;
    BlockIndex block = patch.block;
    auto domain = layout().domain();
    Box domainBox = domain.getBlock(block).box();
    Box patchBox = layout()[index].grow(map.ghost()[0]);
    if (domainBox.containsPoint(pi) && patchBox.containsPoint(pi))
    {
        MBPoint srcPoint(pi, block);

        MBDataPoint data(index, pi, layout());
        sources.push_back(data);
        return true;
    } else if (!domainBox.containsPoint(pi))
    {
        // Examine all patches which border this one across a block boundary
        for (auto dir : Point::Directions())
        {
            if (!layout().isBlockBoundary(index, dir)) { continue; }
            for (auto adjIndex : layout().adjacentIndices(patch, dir))
            {
                BlockIndex adjBlock = layout().block(adjIndex);
                
                Box adjDomainBox = layout().domain().getBlock(adjBlock).box();
                // accessing the bound box using the jacobian is a temporary hack.
                // this data should be queryable in a cleaner fashion
                // jacobian is used since that variable is cell centered (coordinates are node centered)
                Box srcDataBox = map.jacobian().bound(index, dir, adjBlock).adjBox();

                // convert pi -> pj
                // pj is the point associated with the physical coordinate of pi but in the
                // adjacent block's index space
                Point pj = map.convertPoint(pi, block, adjBlock);
                if (srcDataBox.containsPoint(pj))
                {
                    if (!adjDomainBox.containsPoint(pj))
                    {
                        Box dstDataBox = map.map().bound(index, dir, adjBlock).localBox();
                        int _dummy_ = 0;
                    }
                    PROTO_ASSERT(adjDomainBox.containsPoint(pj), "data corruption");
                    // convert pj -> pk
                    // pk is associated with the same logical cellln  as pj but
                    // is in the local block's index space.
                    // NB: In general: pi != pk
                    Point pk = domain.convertPoint(pj, adjBlock, block, PR_CELL);
                    MBDataPoint data(index, pk, adjIndex, pj, layout());
                    sources.push_back(data);
                    return true;
                }
            }
        }
    }
    return false;
}

void MBInterpLayout::removeDomainBoundaryPoints(
    std::set<Point>& footprint,
    Point center,
    MBIndex index) const 
{
    BlockIndex block = m_layout.block(index);
    std::set<Point> tmp(footprint.begin(), footprint.end());
    footprint.clear();
    for (auto p : tmp)
    {
        if (!m_layout.domain().isPointInDomainBoundary(p + center,block)) 
        {
            footprint.insert(p);
        } else {
        }
    }
}


void MBInterpLayout::removeTriplePointPoints(
    std::set<Point>& footprint,
    Point center,
    MBIndex index) const 
{
    BlockIndex block = m_layout.block(index);
    std::set<Point> tmp(footprint.begin(), footprint.end());
    footprint.clear();
    for (auto p : tmp)
    {
        if (!m_layout.domain().isPointInTriplePoint(p + center,block)) 
        {
            footprint.insert(p);
        } else {
        }
    }
}


void MBInterpLayout::addBoundaryFootprint(
    std::set<Point>& footprint,
    Point point,
    MBIndex index) const
{
    auto patch = m_layout.patch(index);
    Box patchBox = m_layout[index];
    std::vector<Point> initialFootprint(footprint.begin(), footprint.end());
    for (auto dir : Point::DirectionsOfCodim(1))
    {
        if (dir * footprintMask() == Point::Zeros()) { continue; }
        if (m_layout.isPatchOnDomainBoundary(patch, dir) ||
            m_layout.isRefinementBoundary(patch, dir))
        {
            int coord = dir.firstNonZeroIndex();
            Box faceBox = patchBox.edge(dir);
            if (point[coord] == faceBox.low()[coord])
            {
                for (auto fi : initialFootprint)
                {
                    footprint.insert(fi - dir);
                }
            }
        }
    }
}


void MBInterpLayout::addTriplePointFootprint(
    std::set<Point>& footprint,
    Point center,
    Point triplePointBoundary,
    MBIndex index) const
{
    BlockIndex block = m_layout.block(index);
    Box patchBox = m_layout[index];

    auto triplePointCircuits = 
        m_layout.domain().graph().getTriplePointCircuits(block, triplePointBoundary);
    PROTO_ASSERT(triplePointCircuits.size() == 1,
        "MBInterpLayout::footprint | Error: The more complex case is not implemented");
    
    // Everything past this point assumes there can only be one triple point circuit

    Point n0 = triplePointCircuits[0].first;
    Point n1 = triplePointCircuits[0].second;
    Point inPlaneMask = (n0 + n1).abs();
    Point outPlaneMask = Point::Ones() - inPlaneMask;
    Point corner = patchBox.edge(n1 + n0).low()*inPlaneMask + center*outPlaneMask;
    int maxDist = (center - corner).abs().max();
    Point s0 = corner + n0*maxDist - center;
    Point s1 = corner + n1*maxDist - center;

    for (auto fi : m_baseFootprint)
    {
        footprint.insert(fi + s0);
        footprint.insert(fi + s1);
    }
}


void MBInterpLayout::addInwardPoints(
        std::set<Point>& footprint,
        MBIndex index) const
{
    std::set<Point> baseFootprint(footprint.begin(), footprint.end());
    Box patchBox = m_layout[index];
    for (auto pi : baseFootprint)
    {
        Point dir = patchBox.whichBoundaryContains(pi);
        footprint.insert(pi - dir);
    }
}


Point MBInterpLayout::footprintMask() const 
{
    Point mask = Point::Ones();
    auto spanBox = footprintSpan();
    for (int dd = 0; dd < DIM; dd++)
    {
        mask[dd] = (spanBox.size(dd) == 1) ? 0 : 1;
    }
    return mask;
}

Box MBInterpLayout::footprintSpan() const 
{
    Box spanBox;
    for (auto fi : m_baseFootprint)
    {
        spanBox &= fi;
    }
    return spanBox;
}