
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
                addPointToBaseFootprint(pi);
            }
        }
        break;
    case 5:
        for (auto pi : Box::Kernel(2))
        {
            if (pi.abs().sum() <= 3)
            {
                addPointToBaseFootprint(pi);
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
        if (fi[dir] == 0) { addPointToBaseFootprint(fi); }
    }
}

template<typename MAP, MemType MEM>
inline std::vector<MBDataPoint>
MBInterpLayout::computeSourcePoints(
    const MBLevelMap<MAP,MEM>& map,
    Point localCenter,
    MBIndex localIndex) const
{
    std::set<MBDataPoint> sourcePoints;
    pr_out() << "\nComputing Source Points | CENTER: " << localCenter << " | PATCH: " << layout().point(localIndex) << " | BLOCK: " << layout().block(localIndex) << std::endl;
    auto trueCenterQuery = map.whichNeighborContains(localIndex, localCenter);
    if (trueCenterQuery.has_value())
    {
        auto [adjIndex, adjCenter] = trueCenterQuery.value();
        Point localTrueCenter = layout().domain().convertPoint(adjCenter, layout().block(adjIndex), layout().block(localIndex), PR_CELL);

        // the following line deals with the case where localTrueCenter has been mapped outside of the hyperplane of the footprint.
        Point projectedLocalTrueCenter = applyIgnoredAxes(localCenter, localTrueCenter);
        pr_out() << "\tComputed true center: " << localTrueCenter << " | projected: " << projectedLocalTrueCenter << std::endl;
        addSourcePointsNear(projectedLocalTrueCenter, localIndex, sourcePoints, map);
        pr_out() << "\tAdded " << sourcePoints.size() << " source points" << std::endl;
    // The else is used for the rare case where a the localCenter doesn't map cleanly to a neighboring patch
    // The best example of this is where 3 mutually symmetric blocks meet
    } else {
        pr_out() << "\tCould not establish a true center. Examining near center points" << std::endl;
        for (auto dir : Point::DirectionsOfCodim(1))
        {
            dir *= footprintMask();
            if (dir == Point::Zeros()) { continue; }
            auto nearCenterQuery = map.whichNeighborContains(localIndex, localCenter + dir);
            if (nearCenterQuery.has_value())
            {
                auto [adjIndex, adjCenter] = nearCenterQuery.value();
                Point localNearCenter = layout().domain().convertPoint(adjCenter, layout().block(adjIndex), layout().block(localIndex), PR_CELL);

                // the following line deals with the case where localTrueCenter has been mapped outside of the hyperplane of the footprint.
                Point projectedLocalNearCenter = applyIgnoredAxes(localCenter, localNearCenter);
                pr_out() << "\tComputed near center: " << localNearCenter << " | projected: " << projectedLocalNearCenter << std::endl;
                int before = sourcePoints.size();
                addSourcePointsNear(projectedLocalNearCenter, localIndex, sourcePoints, map);
                pr_out() << "\tAdded " << sourcePoints.size() - before << " source points" << std::endl;
            }
        }
    }
    // addRedundantBoundarySources(localIndex, sourcePoints, map);
    return std::vector<MBDataPoint>(sourcePoints.begin(), sourcePoints.end());
}

template<typename MAP, MemType MEM>
void MBInterpLayout::addSourcePointsNear(
    Point localTrueCenter,
    MBIndex localIndex,
    std::set<MBDataPoint>& sourcePoints,
    const MBLevelMap<MAP,MEM>& map) const
{
    auto localPatch = layout().patch(localIndex);
    Box localPatchBox = layout()[localIndex];
    std::set<Point> frontierDirs;
    
    // Add all nearby points based on the base footprint
    for (auto fi : baseFootprint())
    {
        Point sourcePoint = localTrueCenter + fi;
        auto adjIndexQuery = map.whichNeighborContains(localIndex, sourcePoint);
        if (adjIndexQuery.has_value())
        {
            auto adjIndex = adjIndexQuery.value().first;
            Point adjPoint = layout().domain().convertPoint(sourcePoint, localPatch.block, layout().block(adjIndex), PR_CELL);
            bool success = addSourcePoint(sourcePoints, sourcePoint, localIndex, adjPoint, adjIndex, map);
            if (!success) { frontierDirs.insert(-(fi.dir())); }
        } else {
            // if unable to add this point, store the opposite direction for extending the footprint away from the "hole"
            frontierDirs.insert(-(fi.dir()));
        }
    }

    // if necessary, process an additional frontier on the edge of the footprint
    if (frontierDirs.size() > 0)
    {
        auto tmp = baseFootprint();
        std::set<Point> footprint(tmp.begin(), tmp.end());
        std::set<Point> frontier;
        for (auto fi : footprint)
        {
            for (auto dir : Point::Directions())
            {
                Point fj = fi + dir;
                if (frontierDirs.count(dir) > 0 && footprint.count(fj) <= 0)
                {
                    frontier.insert(fi + dir);
                }
            }
        }

        // this loop is the same as the footprint loop except that points which don't map properly are simply ignored
        for (auto fi : frontier)
        {
            Point sourcePoint = localTrueCenter + fi;
            auto adjIndexQuery = map.whichNeighborContains(localIndex, sourcePoint);
            if (adjIndexQuery.has_value())
            {
                auto adjIndex = adjIndexQuery.value().first;
                Point adjPoint = layout().domain().convertPoint(sourcePoint, localPatch.block, layout().block(adjIndex), PR_CELL);
                addSourcePoint(sourcePoints, sourcePoint, localIndex, adjPoint, adjIndex, map);
            }
        }
    }
}

template<typename MAP, MemType MEM>
bool MBInterpLayout::addSourcePoint(
    std::set<MBDataPoint>& sourcePoints,
    Point localPoint,
    MBIndex localIndex,
    Point sourcePoint,
    MBIndex sourceIndex,
    const MBLevelMap<MAP,MEM>& map) const
{
    Box localPatchBox = layout()[localIndex];
    Box localDomainBox = layout().domain().getBlock(layout().block(localIndex)).box();
    if (layout().block(localIndex) == layout().block(sourceIndex))
    {
        if ((localPatchBox.grow(map.ghost()[0]) & localDomainBox).containsPoint(localPoint))
        {
            MBDataPoint point(localIndex, localPoint, layout());
            auto [_, success] = sourcePoints.insert(point);
            return success;
        }
        pr_out() << "FAILURE\n";
        // otherwise point is in an adjacent patch in the same block but NOT IN THE GHOST REGION of this patch and we ignore it
    } else {
        Point boundaryDir = layout().connectivity(localIndex, sourceIndex);
        Point boundGhost = map.ghost()[boundaryDir.codim()];
        Box boundaryBox = localPatchBox.adjacent(boundaryDir * boundGhost);
        if ( boundaryDir != Point::Zeros() && boundaryBox.containsPoint(localPoint))
        {
            MBDataPoint point(localIndex, localPoint, sourceIndex, sourcePoint, layout());
            auto [_, success] = sourcePoints.insert(point);
            return success;
        }
    }
    return false;
}

template<typename MAP, MemType MEM>
void MBInterpLayout::addRedundantBoundarySources(
    MBIndex localIndex,
    std::set<MBDataPoint>& sourcePoints,
    const MBLevelMap<MAP, MEM>& map) const
{
    auto localPatch = layout().patch(localIndex);
    std::set<Point> localSourcePoints;
    for (auto& si : sourcePoints) { localSourcePoints.insert(si.dstPoint()); }
    Box localPatchBox = layout()[localIndex];
    for (auto boundaryDir : Point::Directions())
    {
        if (boundaryDir.codim() < 2) { continue; }
        Box nearBoundaryThreshold = localPatchBox.edge(boundaryDir).extrude(boundaryDir);
        for (auto p0 : localSourcePoints)
        {
            if (nearBoundaryThreshold.containsPoint(p0))
            {
                Box nearPointThreshold(p0, p0);
                nearPointThreshold = nearPointThreshold.extrude(boundaryDir * footprintMask());
                for (auto dir : boundaryDir.subset())
                {
                    for (auto adjIndex : layout().adjacentIndices(localPatch, dir))
                    {
                        auto adjBlock = layout().block(adjIndex);
                        Box nearPointThresholdAdj = layout().domain().convertBox(nearPointThreshold, localPatch.block, adjBlock);
                        nearPointThresholdAdj &= layout()[adjIndex];
                        for (auto srcPointAdj : nearPointThresholdAdj)
                        {
                            Point srcPointLocal = layout().domain().convertPoint(srcPointAdj, adjBlock, localPatch.block, PR_CELL);
                            bool success = addSourcePoint(sourcePoints, srcPointLocal, localIndex, srcPointAdj, adjIndex, map);
                        }
                    }
                }
            }
        }
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

Point MBInterpLayout::applyIgnoredAxes(Point basePoint, Point sourcePoint) const
{
    Point diff = sourcePoint - basePoint;
    Point mask = Point::Ones() - footprintMask();
    diff = diff * mask;
    return sourcePoint - diff;
}