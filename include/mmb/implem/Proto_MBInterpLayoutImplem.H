
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

    auto trueCenterQuery = map.whichNeighborContains(localIndex, localCenter);
    if (trueCenterQuery.has_value())
    {
        auto [adjIndex, adjCenter] = trueCenterQuery.value();
        Point localTrueCenter = layout().domain().convertPoint(adjCenter, layout().block(adjIndex), layout().block(localIndex), PR_CELL);
        addSourcePointsNear(localTrueCenter, localIndex, sourcePoints, map);
    // The else is used for the rare case where a the localCenter doesn't map cleanly to a neighboring patch
    // The best example of this is where 3 mutually symmetric blocks meet
    } else {
        for (auto dir : Point::DirectionsOfCodim(1))
        {
            auto nearCenterQuery = map.whichNeighborContains(localIndex, localCenter + dir);
            if (nearCenterQuery.has_value())
            {
                auto [adjIndex, adjCenter] = nearCenterQuery.value();
                Point localNearCenter = layout().domain().convertPoint(adjCenter, layout().block(adjIndex), layout().block(localIndex), PR_CELL);
                addSourcePointsNear(localNearCenter, localIndex, sourcePoints, map);
            }
        }
    }
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
        pr_out() << "\t\t\tTrying to add source from local patch...";
        if ((localPatchBox.grow(map.ghost()[0]) & localDomainBox).containsPoint(localPoint))
        {
            MBDataPoint point(localIndex, localPoint, layout());
            auto [_, success] = sourcePoints.insert(point);

            if (success) { pr_out() << "SUCCESS\n"; } else { pr_out() << "FAILURE\n"; }
            return success;
        }
        pr_out() << "FAILURE\n";
        // otherwise point is in an adjacent patch in the same block but NOT IN THE GHOST REGION of this patch and we ignore it
    } else {
        Point boundaryDir = layout().connectivity(localIndex, sourceIndex);
        Point boundGhost = map.ghost()[boundaryDir.codim()];
        Box boundaryBox = localPatchBox.adjacent(boundaryDir * boundGhost);
        pr_out() << "\t\t\tTrying to add source from adjacent block...";
        if ( boundaryDir != Point::Zeros() && boundaryBox.containsPoint(localPoint))
        {
            MBDataPoint point(localIndex, localPoint, sourceIndex, sourcePoint, layout());
            auto [_, success] = sourcePoints.insert(point);
            if (success) { pr_out() << "SUCCESS\n"; } else { pr_out() << "FAILURE\n"; }
            return success;
        }
        pr_out() << "FAILURE\n";
    }
    return false;
}

// void MBInterpLayout::removeDomainBoundaryPoints(
//     std::set<Point>& footprint,
//     Point center,
//     MBIndex index) const 
// {
//     BlockIndex block = m_layout.block(index);
//     std::set<Point> tmp(footprint.begin(), footprint.end());
//     footprint.clear();
//     for (auto p : tmp)
//     {
//         if (!m_layout.domain().isPointInDomainBoundary(p + center,block)) 
//         {
//             footprint.insert(p);
//         } else {
//         }
//     }
// }


// void MBInterpLayout::removeTriplePointPoints(
//     std::set<Point>& footprint,
//     Point center,
//     MBIndex index) const 
// {
//     BlockIndex block = m_layout.block(index);
//     std::set<Point> tmp(footprint.begin(), footprint.end());
//     footprint.clear();
//     for (auto p : tmp)
//     {
//         if (!m_layout.domain().isPointInTriplePoint(p + center,block)) 
//         {
//             footprint.insert(p);
//         } else {
//         }
//     }
// }


// void MBInterpLayout::addBoundaryFootprint(
//     std::set<Point>& footprint,
//     Point point,
//     MBIndex index) const
// {
//     auto patch = m_layout.patch(index);
//     Box patchBox = m_layout[index];
//     std::vector<Point> initialFootprint(footprint.begin(), footprint.end());
//     for (auto dir : Point::DirectionsOfCodim(1))
//     {
//         if (dir * footprintMask() == Point::Zeros()) { continue; }
//         if (m_layout.isPatchOnDomainBoundary(patch, dir) ||
//             m_layout.isRefinementBoundary(patch, dir))
//         {
//             int coord = dir.firstNonZeroIndex();
//             Box faceBox = patchBox.edge(dir);
//             if (point[coord] == faceBox.low()[coord])
//             {
//                 for (auto fi : initialFootprint)
//                 {
//                     footprint.insert(fi - dir);
//                 }
//             }
//         }
//     }
// }


// void MBInterpLayout::addTriplePointFootprint(
//     std::set<Point>& footprint,
//     Point center,
//     Point triplePointBoundary,
//     MBIndex index) const
// {
//     BlockIndex block = m_layout.block(index);
//     Box patchBox = m_layout[index];

//     auto triplePointCircuits = 
//         m_layout.domain().graph().getTriplePointCircuits(block, triplePointBoundary);
//     PROTO_ASSERT(triplePointCircuits.size() == 1,
//         "MBInterpLayout::footprint | Error: The more complex case is not implemented");
    
//     // Everything past this point assumes there can only be one triple point circuit

//     Point n0 = triplePointCircuits[0].first;
//     Point n1 = triplePointCircuits[0].second;
//     Point inPlaneMask = (n0 + n1).abs();
//     Point outPlaneMask = Point::Ones() - inPlaneMask;
//     Point corner = patchBox.edge(n1 + n0).low()*inPlaneMask + center*outPlaneMask;
//     int maxDist = (center - corner).abs().max();
//     Point s0 = corner + n0*maxDist - center;
//     Point s1 = corner + n1*maxDist - center;

//     for (auto fi : m_baseFootprint)
//     {
//         footprint.insert(fi + s0);
//         footprint.insert(fi + s1);
//     }
// }


// void MBInterpLayout::addInwardPoints(
//         std::set<Point>& footprint,
//         MBIndex index) const
// {
//     std::set<Point> baseFootprint(footprint.begin(), footprint.end());
//     Box patchBox = m_layout[index];
//     for (auto pi : baseFootprint)
//     {
//         Point dir = patchBox.whichBoundaryContains(pi);
//         footprint.insert(pi - dir);
//     }
// }


// Point MBInterpLayout::footprintMask() const 
// {
//     Point mask = Point::Ones();
//     auto spanBox = footprintSpan();
//     for (int dd = 0; dd < DIM; dd++)
//     {
//         mask[dd] = (spanBox.size(dd) == 1) ? 0 : 1;
//     }
//     return mask;
// }

// Box MBInterpLayout::footprintSpan() const 
// {
//     Box spanBox;
//     for (auto fi : m_baseFootprint)
//     {
//         spanBox &= fi;
//     }
//     return spanBox;
// }