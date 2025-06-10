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

        inline void addInwardPoints(
            std::set<Point>& footprint,
            MBIndex index) const;

        MBDisjointBoxLayout m_layout;
        std::set<Point> m_baseFootprint;
        Array<bool, DIM> m_copyAxes;
    };
    #include "implem/Proto_MBInterpLayoutImplem.H"

} //end namespace Proto
#endif // end include guard
