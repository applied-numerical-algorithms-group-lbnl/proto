
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
    
    //find which block the center corresponds to
        // need to try all blocks with ANY boundary with this PATCH
    // find center cell in adjacent block
        // SHOULD be unique
    auto localPatch = layout().patch(localIndex);
    Box localPatchBox = layout()[localIndex];
    Box localDomainBox = layout().domain().getBlock(localPatch.block).box();
    Point adjCenter;
    MBIndex adjIndex;
    PROTO_ASSERT(!localDomainBox.containsPoint(localCenter), "Trying to interpolate to an interior point");
    // std::cout << "Trying to find the block that contains localCenter " << localCenter << std::endl;
    for (auto neighbor : layout().adjacentIndices(localPatch))
    {
        Box neighborBox = layout()[neighbor];
        if (localCenter == Point(16,16))
        {
            int _dummy = 0;
        }
        adjCenter = map.convertPoint(localCenter, localPatch.block, layout().block(neighbor));
        // std::cout << "Computed adjCenter: " << adjCenter << " | neighborBox: " << neighborBox << std::endl;
        if (neighborBox.containsPoint(adjCenter))
        {
            // std::cout << "Found it. Maps to adjCenter = " << adjCenter << " in block " << layout().block(neighbor) << std::endl;
            adjIndex = neighbor;
            break;
        }
    }
    PROTO_ASSERT(layout()[adjIndex].containsPoint(adjCenter),
        "Center point does not map into an adjacent block. You might be trying to access data in a DOMAIN boundary.");
    MBPoint adjPatch = layout().patch(adjIndex);
    Box adjPatchBox = layout()[adjIndex];

    pr_out() << "\nComputing sources for POINT: " << localCenter << " | BLOCK: " << localPatch.block << " | PATCH: " << localPatch.point << std::endl;
    pr_out() << "\tCenter is mapped to POINT: " << adjCenter << " | BLOCK: " << adjPatch.block << " | PATCH:" << adjPatch.point << std::endl;
    std::set<MBDataPoint> sourcePoints;
    std::set<Point> extendDirs;

    Point localTrueCenter = layout().domain().convertPoint(adjCenter, adjPatch.block, localPatch.block, PR_CELL);
    for (auto fi : baseFootprint())
    {
        Point sourcePoint = localTrueCenter + fi;
        bool success = false;
        pr_out() << "\t\tTrying to add source point " << sourcePoint << " from shift " << fi;
        if (localPatchBox.containsPoint(sourcePoint))
        {
            success = addSourcePoint(sourcePoints, sourcePoint, localIndex, sourcePoint, localIndex, map);
            if (success)
            {
                pr_out() << "\t\tSUCCESS (found in local patch)" << std::endl;
            }
        } else {
            for (auto adjIndex_ : layout().adjacentIndices(localPatch))
            {
                MBPoint adjPatch_ = layout().patch(adjIndex_);
                Point adjSourcePoint = map.convertPoint(sourcePoint, localPatch.block, adjPatch_.block);
                Box adjPatchBox = layout()[adjIndex_];
                if (adjPatchBox.containsPoint(adjSourcePoint))
                {
                    success = addSourcePoint(sourcePoints, sourcePoint, localIndex, adjSourcePoint, adjIndex_, map);
                    if (success)
                    {
                        pr_out() << "\t\tSUCCESS (found in adjacent patch)" << std::endl;
                    }
                    break;
                }
            }
        }
        if (!success)
        {
            pr_out() << "\t\tFAILURE (point is in a domain boundary)" << std::endl;
            extendDirs.insert(-(fi.dir()));
        }
    }

    pr_out() << "Extending footprint in directions: " << std::endl;
    for (auto dir : extendDirs) { pr_out() << dir << ", "; } pr_out() << std::endl;

    if (extendDirs.size() > 0)
    {
        auto tmp = baseFootprint();
        std::set<Point> footprint(tmp.begin(), tmp.end());
        for (auto fi : baseFootprint())
        {
            for (auto dir : Point::Directions())
            {
                if (extendDirs.count(dir) > 0)
                {
                    footprint.insert(fi + dir);
                }
            }
        }

        for (auto fi : footprint)
        {
            Point sourcePoint = localTrueCenter + fi;
            bool success = false;
            pr_out() << "\t\tTrying to add source point " << sourcePoint << " from shift " << fi;
            if (localPatchBox.containsPoint(sourcePoint))
            {
                success = addSourcePoint(sourcePoints, sourcePoint, localIndex, sourcePoint, localIndex, map);
                if (success)
                {
                    pr_out() << "\t\tSUCCESS (found in local patch)" << std::endl;
                }
            } else {
                for (auto adjIndex_ : layout().adjacentIndices(localPatch))
                {
                    MBPoint adjPatch_ = layout().patch(adjIndex_);
                    Point adjSourcePoint = map.convertPoint(sourcePoint, localPatch.block, adjPatch_.block);
                    Box adjPatchBox = layout()[adjIndex_];
                    if (adjPatchBox.containsPoint(adjSourcePoint))
                    {
                        success = addSourcePoint(sourcePoints, sourcePoint, localIndex, adjSourcePoint, adjIndex_, map);
                        if (success)
                        {
                            pr_out() << "\t\tSUCCESS (found in adjacent patch)" << std::endl;
                        }
                        break;
                    }
                }
            }
            if (!success)
            {
                pr_out() << "\t\tFAILURE (point is duplicate or in a domain boundary)" << std::endl;
            }
        }
    }

    // for (auto fi : footprint)
    // {
    //     Point fj = layout().domain().convertPoint(fi, localPatch.block, adjPatch.block, PR_CELL);
    //     Point adjSourcePoint = adjCenter + fj;
    //     pr_out() << "\tTrying to add POINT: " << adjSourcePoint << "\n";
    //     if (adjPatchBox.containsPoint(adjSourcePoint))
    //     {
    //         Point localSourcePoint = map.convertPoint(adjSourcePoint, adjPatch.block, localPatch.block);
    //         bool success = addSourcePoint(sourcePoints, localSourcePoint, localIndex, adjSourcePoint, adjIndex, map);
    //         //if (!success){ extendDirs.insert(-fi.dir()); }
    //         pr_out() << "\t\tPoint found | Maps to " << localSourcePoint << " | block: " << adjPatch.block << " | patch: " << adjPatch.point;
    //         if (success)
    //         {
    //             pr_out() << "\t\tSUCCESS" << std::endl;
    //         } else {
    //             pr_out() << "\t\tFAILURE" << std::endl;
    //         }
    //     } else {
    //         for (auto adjAdjIndex : layout().adjacentIndices(adjPatch))
    //         {
    //             MBPoint adjAdjPatch = layout().patch(adjAdjIndex);
    //             Point adjAdjSourcePoint = map.convertPoint(adjSourcePoint, adjPatch.block, adjAdjPatch.block);
    //             Box adjAdjPatchBox = layout()[adjAdjIndex];
    //             if (adjAdjPatchBox.containsPoint(adjAdjSourcePoint))
    //             {
    //                 Point localSourcePoint = map.convertPoint(adjAdjSourcePoint, adjAdjPatch.block, localPatch.block);
    //                 bool success = addSourcePoint(sourcePoints, localSourcePoint, localIndex, adjAdjSourcePoint, adjAdjIndex, map);
    //                 //if (!success) { extendDirs.insert(-fi.dir()); }
    //                 pr_out() << "\t\tPoint found | Maps to " << localSourcePoint << " | block: " << adjAdjPatch.block << " | patch: " << adjAdjPatch.point;
    //                 if (success)
    //                 {
    //                     pr_out() << "\t\tSUCCESS" << std::endl;
    //                 } else {
    //                     pr_out() << "\t\tFAILURE" << std::endl;
    //                 }
    //                 break;
    //             }
    //         }
    //     }
    // }


    return std::vector<MBDataPoint>(sourcePoints.begin(), sourcePoints.end());
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
        if ((localDomainBox.grow(map.ghost()[0]) & localPatchBox).containsPoint(localPoint))
        {
            MBDataPoint point(localIndex, localPoint, layout());
            auto [_, success] = sourcePoints.insert(point);
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
            return success;
        }
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