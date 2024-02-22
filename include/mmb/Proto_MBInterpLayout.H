#pragma once
#ifndef _PROTO_MB_INTERP_LAYOUT_
#define _PROTO_MB_INTERP_LAYOUT_

#include "Proto.H"

namespace Proto
{
    class MBInterpLayout
    {
        public:
        
        inline MBInterpLayout(){ m_footprint.clear(); m_copyAxes.fill(false); m_mirrorAxes.fill(false); }
        inline bool addPoint(Point a_point) { m_footprint.insert(a_point); }
        inline void setCopy(int a_dir, bool a_flag) { m_copyAxes[a_dir] = a_flag; if (a_flag) {m_mirrorAxes[a_dir] = false; }}
        inline void setMirror(int a_dir, bool a_flag) { m_mirrorAxes[a_dir] = a_flag; if (a_flag) {m_copyAxes[a_dir] = false; }}

        inline bool doMirror(int a_dir) const { return m_mirrorAxes[a_dir]; }
        inline bool doCopy(int a_dir) const { return m_copyAxes[a_dir]; }
        inline std::vector<Point> footprint() const { return std::vector<Point>(m_footprint.begin(), m_footprint.end()); }
        
        private:
        std::set<Point> m_footprint;
        Array<bool, DIM> m_copyAxes;
        Array<bool, DIM> m_mirrorAxes;
    };

//#include "implem/Proto_MBInterpLayoutImplem.H"
} //end namespace Proto
#endif // end include guard
