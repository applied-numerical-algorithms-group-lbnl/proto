#pragma once
#ifndef _PROTO_MB_INTERP_LAYOUT_
#define _PROTO_MB_INTERP_LAYOUT_

#include "Proto.H"

namespace Proto
{
    class MBInterpLayout
    {
        public:
        inline MBInterpLayout(){}
        inline MBInterpLayout(const MBDisjointBoxLayout& layout, int order = 4);
        inline void addPoint(Point a_point) { m_baseFootprint.insert(a_point); }
        inline void ignoreAxis(int dir);
        inline void setCopy(int a_dir, bool a_flag) { m_copyAxes[a_dir] = a_flag; }

        inline bool isCopied(int a_dir) const { return m_copyAxes[a_dir]; }
        inline std::set<Point> baseFootprint() const { return m_baseFootprint; }
        inline std::set<Point> footprint(
            Point point,
            Point boundaryDir,
            MBIndex index) const;
        inline Point footprintMask() const;
        inline Box footprintSpan() const;
        inline const MBDisjointBoxLayout& layout() const { return m_layout; }

        private:

        inline void addTriplePointFootprint(
            std::set<Point>& footprint,
            Point center,
            Point triplePointBoundary,
            MBIndex index) const;
        
        inline void addBoundaryFootprint(
            std::set<Point>& footprint,
            Point center,
            MBIndex index) const;

        inline void removeDomainBoundaryPoints(
            std::set<Point>& footprint,
            Point center,
            MBIndex index) const;

        inline void removeTriplePointPoints(
            std::set<Point>& footprint,
            Point center,
            MBIndex index) const;

        MBDisjointBoxLayout m_layout;
        std::set<Point> m_baseFootprint;
        Array<bool, DIM> m_copyAxes;
    };

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
        Point boundaryDir,
        MBIndex index) const
    {
        PROTO_ASSERT(boundaryDir.abs().max() <= 1,
            "MBInterpLayout::footprint | Error: boundaryDir is not a valid direction");
        auto patchDomain = m_layout.patchDomain();
        BlockIndex block = m_layout.block(index);
        PatchID patch = m_layout.point(index);
        PatchID adjPatch = patch + boundaryDir;
        Box patchBox = m_layout[index];
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

} //end namespace Proto
#endif // end include guard
