#pragma once
#ifndef _PROTO_MB_INTERP_LAYOUT_
#define _PROTO_MB_INTERP_LAYOUT_

#include "Proto.H"

namespace Proto
{
    class MBInterpLayout
    {
        public:
        
        inline MBInterpLayout(MBDisjointBoxLayout& layout);
        inline bool addPoint(Point a_point) { m_baseFootprint.insert(a_point); }
        inline void setCopy(int a_dir, bool a_flag) { m_copyAxes[a_dir] = a_flag; }

        inline bool isCopied(int a_dir) const { return m_copyAxes[a_dir]; }
        inline std::set<Point> baseFootprint() const { return m_baseFootprint; }
        inline std::set<Point> footprint(
            Point point,
            Point boundaryDir,
            MBIndex index) const;
        inline Point footprintMask() const;
        inline Box footprintSpan() const;

        private:

        inline void addTriplePointFootprint(
            std::set<Point>& footprint,
            Point center,
            Point triplePointBoundary,
            MBIndex index) const;
        
        inline void addDomainBoundaryFootprint(
            std::set<Point>& footprint,
            Point point,
            MBIndex index) const;

        MBDisjointBoxLayout m_layout;
        std::set<Point> m_baseFootprint;
        Array<bool, DIM> m_copyAxes;
    };

    MBInterpLayout::MBInterpLayout(MBDisjointBoxLayout& layout)
    {
        m_baseFootprint.clear();
        m_copyAxes.fill(false);
        m_layout = layout;
    }

    std::set<Point> MBInterpLayout::footprint(
        Point point,
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
        addDomainBoundaryFootprint(footprint, point, index);
        if (patchDomain.inTriplePoint(adjPatch, block))
        {
            addTriplePointFootprint(footprint, point, boundaryDir, index);
        } else if (patchDomain.inBlockBoundary(adjPatch, block))
        {
            Point orthDir = Point::Ones() - adjPatch.abs();
            std::vector<Point> orthDirs;
            for (auto orth : orthDir.split())
            {
                orthDirs.push_back(orth);
                orthDirs.push_back(-orth);
            }
            for (auto orth : orthDirs)
            {
                if (patchDomain.inTriplePoint(adjPatch + orth, block))
                {
                    addTriplePointFootprint(footprint, point, boundaryDir, index);
                }
            }
        }
        return footprint;
    }

    void MBInterpLayout::addDomainBoundaryFootprint(
        std::set<Point>& footprint,
        Point point,
        MBIndex index) const
    {
        PatchID patch = m_layout.point(index);
        BlockIndex block = m_layout.block(index);
        Box patchBox = m_layout[index];
        std::vector<Point> initialFootprint(footprint.begin(), footprint.end());
        for (auto dir : Point::DirectionsOfCodim(1))
        {
            if (dir * footprintMask() == Point::Zeros()) { continue; }
            if (m_layout.onDomainBoundary(patch, block, dir))
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
        Point inPlane = n0 + n1; // for DIM=2, this is the same as boundaryDir
        Point outPlane = (Point::Ones() - inPlane.abs())*triplePointBoundary;
        Point corner = patchBox.edge(inPlane).low()*inPlane;
        Point posFromCorner = (center - corner)*inPlane;
        int maxDistFromCorner = posFromCorner.abs().max();
        Point shift0 = corner + n0*maxDistFromCorner + center*outPlane - center;
        Point shift1 = corner + n1*maxDistFromCorner + center*outPlane - center;

        for (auto fi : m_baseFootprint)
        {
            Point f0 = fi + shift0;
            Point f1 = fi + shift1;
            footprint.insert(f0);
            footprint.insert(f1);
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
