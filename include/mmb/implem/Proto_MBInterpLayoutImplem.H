
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
    MBIndex localIndex,
    int minSources) const
{
    std::set<MBDataPoint> sourcePoints;
    pr_out() << "\nComputing Source Points | CENTER: " << localCenter << " | PATCH: " << layout().point(localIndex) << " | BLOCK: " << layout().block(localIndex) << std::endl;
    auto trueCenterQuery = map.whichNeighborContains(localIndex, localCenter);
    Point trueCenter = Point::Zeros();
    if (trueCenterQuery.has_value())
    {
        auto [adjIndex, adjCenter] = trueCenterQuery.value();

        // the following line deals with the case where localTrueCenter has been mapped outside of the hyperplane of the footprint.
        Point projectedAdjCenter = applyIgnoredAxes(localCenter, adjCenter);
        addSourcePointsNear(projectedAdjCenter, adjIndex, localIndex, sourcePoints, map);
        Point localTrueCenter = layout().domain().convertPoint(projectedAdjCenter, layout().block(adjIndex), layout().block(localIndex), PR_CELL);
        addSourcePointsNear(localTrueCenter, localIndex, localIndex, sourcePoints, map);
        trueCenter = localTrueCenter;
    // The else is used for the rare case where a the localCenter doesn't map cleanly to a neighboring patch
    // The best example of this is where an odd number of mutually symmetric blocks meet
    } else {
        pr_out() << "\tCould not establish a true center. Examining near center points" << std::endl;
        for (auto dir : Point::DirectionsOfCodim(1))
        {
            if (dir * footprintMask() == Point::Zeros()) { continue; }
            auto nearCenterQuery = map.whichNeighborContains(localIndex, localCenter + dir);
            if (nearCenterQuery.has_value())
            {
                auto [adjIndex, adjCenter] = nearCenterQuery.value();

                // the following line deals with the case where localTrueCenter has been mapped outside of the hyperplane of the footprint.
                Point projectedAdjNearCenter = applyIgnoredAxes(localCenter, adjCenter);
                addSourcePointsNear(projectedAdjNearCenter, adjIndex, localIndex, sourcePoints, map);
                Point localNearCenter = layout().domain().convertPoint(projectedAdjNearCenter, layout().block(adjIndex), layout().block(localIndex), PR_CELL);
                addSourcePointsNear(localNearCenter, localIndex, localIndex, sourcePoints, map);
                trueCenter = (trueCenter == Point::Zeros()) ? localNearCenter : (trueCenter + localNearCenter) / Point::Ones(2);
            }
        }
    }
    // addRedundantBoundarySources is a subroutine which is intended to ensure that 
    // all blocks which share a codimension 2+ boundary contribute at least 1 source point
    // for interpolations which are close enough to said boundary. This theoretically
    // can improve solver stability, however there has not been any need for it as of yet
    // so it is disabled for now - Chris Bozhart - 8/26/2025
    // addRedundantBoundarySources(localIndex, sourcePoints, map);

    // remove excessive source points to prevent loss of accuracy
    // I have found that it makes sense to limit the sources to 2*minSources
    // but it's possible that this value is not optimal. - Chris Bozhart - 8/26/2025
    removeSourcePointsNear(trueCenter, 2*minSources, sourcePoints);
    return std::vector<MBDataPoint>(sourcePoints.begin(), sourcePoints.end());
}

template<typename MAP, MemType MEM>
void MBInterpLayout::addSourcePointsNear(
    Point centerPoint,
    MBIndex searchIndex,
    MBIndex localIndex,
    std::set<MBDataPoint>& sourcePoints,
    const MBLevelMap<MAP,MEM>& map) const
{
    auto localPatch = layout().patch(localIndex);
    auto searchPatch = layout().patch(searchIndex);
    std::set<Point> frontierDirs;

    for (auto fi : baseFootprint())
    {
        Point sourcePoint = centerPoint + fi;
        auto adjIndexQuery = map.whichNeighborContains(searchIndex, sourcePoint);
        if (adjIndexQuery.has_value())
        {
            auto [adjIndex, adjPoint] = adjIndexQuery.value();
            
            if (layout().connectivity(localIndex, adjIndex) != Point::Zeros() || 
                adjIndex == localIndex)
            {
                // sourcePoint maps to localPatch or a patch which shares a boundary with it
                Point sourcePointLocal = layout().domain().convertPoint(adjPoint, layout().block(adjIndex), localPatch.block, PR_CELL);
                bool success = addSourcePoint(sourcePoints, sourcePointLocal, localIndex, adjPoint, adjIndex, map);
                if (!success) { frontierDirs.insert(-(fi.dir())); }

            } else
            {
                // sourcePoint maps to a patch other than localPatch which isn't connected to localPatch
                frontierDirs.insert(-(fi.dir())); 
            }

        } else {
            // sourcePoint doesn't map to a patch (e.g outside of boundary)
            frontierDirs.insert(-fi.dir());
        }
    }

    // if necessary, process an additional frontier on the edge of the footprint
    // for difficult maps, it may be necessary to process two frontiers, but for
    // general use, one is sufficient
    int numFrontiers = 2;
    auto tmp = baseFootprint();
    std::set<Point> footprint(tmp.begin(), tmp.end());
    if (frontierDirs.size() > 0)
    {
        for (int ii = 0; ii < numFrontiers; ii++)
        {
            std::set<Point> frontier;
            for (auto fi : footprint)
            {
                for (auto dir : Point::Directions())
                {
                    if (dir * footprintMask() == Point::Zeros()) { continue; }
                    Point fj = fi + dir;
                    if (frontierDirs.count(dir) > 0 && footprint.count(fj) <= 0)
                    {
                        frontier.insert(fj);
                    }
                }
            }

            // this loop is the same as the footprint loop except that points which don't map properly are simply ignored
            for (auto fi : frontier)
            {
                Point sourcePoint = centerPoint + fi;
                auto adjIndexQuery = map.whichNeighborContains(searchIndex, sourcePoint);
                if (adjIndexQuery.has_value())
                {
                    auto [adjIndex, adjPoint] = adjIndexQuery.value();
                    
                    if (layout().connectivity(localIndex, adjIndex) != Point::Zeros() || adjIndex == localIndex)
                    {
                        // sourcePoint maps to localPatch or a patch which shares a boundary with it
                        Point sourcePointLocal = layout().domain().convertPoint(adjPoint, layout().block(adjIndex), localPatch.block, PR_CELL);
                        bool success = addSourcePoint(sourcePoints, sourcePointLocal, localIndex, adjPoint, adjIndex, map);
                    }
                }
                footprint.insert(fi);
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
    // pr_out() << "trying to add source point | local point: " << localPoint << " | block: " << layout().block(localIndex);
    // pr_out() << " | adj point: " << sourcePoint << " | adj block: " << layout().block(sourceIndex);
    if (layout().block(localIndex) == layout().block(sourceIndex))
    {
        if ((localPatchBox.grow(map.ghost()[0]) & localDomainBox).containsPoint(localPoint))
        {
            MBDataPoint point(localIndex, localPoint, layout());
            auto [_, success] = sourcePoints.insert(point);
            // if (success) { pr_out() << "\tSUCCESS\nAdded point from local patch"<< std::endl;  } else { pr_out() << "\tFAILURE\nPoint is a duplicate"<< std::endl;  }
            return success;
        }
        // otherwise point is in an adjacent patch in the same block but NOT IN THE GHOST REGION of this patch and we ignore it
    } else {
        Point boundaryDir = layout().connectivity(localIndex, sourceIndex);
        Point boundGhost = map.ghost()[boundaryDir.codim()];
        Box boundaryBox = localPatchBox.adjacent(boundaryDir * boundGhost);
        if ( boundaryDir != Point::Zeros() && boundaryBox.containsPoint(localPoint))
        {
            MBDataPoint point(localIndex, localPoint, sourceIndex, sourcePoint, layout());
            auto [_, success] = sourcePoints.insert(point);
            // if (success) { pr_out() << "\tSUCCESS\nAdded point from adjacent patch" << std::endl; } else { pr_out() << "\tFAILURE\nPoint is a duplicate"<< std::endl;  }
            return success;
        }
    }
    // pr_out() << "\tFAILURE\nPoint is not accessible from the local patch"<< std::endl; 
    return false;
}

void MBInterpLayout::removeSourcePointsNear(
    Point trueCenterPoint,
    int maxSources,
    std::set<MBDataPoint>& sourcePoints) const
{
    if (sourcePoints.size() <= maxSources) { return;}
    std::vector<MBDataPoint> sources(sourcePoints.begin(), sourcePoints.end());
    std::sort(sources.begin(), sources.end(),
        [trueCenterPoint](MBDataPoint& a, MBDataPoint& b)
        {
            Point p = a.dstPoint() - trueCenterPoint;
            Point q = b.dstPoint() - trueCenterPoint;
            int dp = 0; int dq = 0;
            for (int dir = 0; dir < DIM; dir++)
            {
                dp += p[dir]*p[dir];
                dq += q[dir]*q[dir];
            }
            return dp < dq;
        });
    
    sourcePoints.clear();
    for (int ii = 0; ii < maxSources; ii++)
    {
        sourcePoints.insert(sources[ii]);
    }
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