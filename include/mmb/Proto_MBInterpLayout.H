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

        /// @brief Modifies the base footprint so that the designated axis is ignored
        /// @param dir The axis to be ignored
        inline void ignoreAxis(int dir);

        /// @brief Allows stencils produced with this layout to be copied along the designated axis (optimization)
        /// @param a_dir 
        /// @param a_flag 
        inline void setCopy(int a_dir, bool a_flag) { m_copyAxes[a_dir] = a_flag; }

        inline bool isCopied(int a_dir) const { return m_copyAxes[a_dir]; }

        /// @brief Returns the set of shifts that would be used in the bulk of a regular grid
        /// @return 
        inline std::set<Point> baseFootprint() const { return m_baseFootprint; }

        /// @brief Returns the modified base footprint, adding extra points for dealing with boundaries
        /// @param point The "center" of the stencil. This is the point being interpolated to
        /// @param index The index of the patch which contains point (in its ghost region)
        /// @return 
        inline std::set<Point> footprint(
            Point point,
            MBIndex index) const;
        
        /// @brief Converts a set of shifts (footprint) into a set of MBDataPoint used to access data for interpolation
        /// @param footprint The set of shifts returned from this->footprint above
        /// @param map Used to convert local index space points into points in adjacent blocks as needed
        /// @param center The point being interpolated to
        /// @param index The index of the patch containing center
        /// @return 
        template<typename MAP, MemType MEM>
        inline std::vector<MBDataPoint> computeSourcePoints(
            std::set<Point> footprint,
            const MBLevelMap<MAP, MEM>& map,
            Point center, MBIndex index) const;



        /// @brief Return contains 1 for axes participating in interpolation and 0 for axes not participating
        /// @return 
        inline Point footprintMask() const;

        /// @brief Returns the span of the base footprint (for 4th order this should be Box::Kernel(2) if no axes are ignored)
        /// @return 
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

        // Currently unused
        inline void addInwardPoints(
            std::set<Point>& footprint,
            MBIndex index) const;

        inline void augmentFootprint(
            Point shift,
            std::set<Point>& augmentedFootprint,
            const std::set<Point>& footprint) const;

        template<typename MAP, MemType MEM>
        inline bool addSourceFromShiftedPoint(
            Point shiftedPoint,
            std::set<MBDataPoint>& sources,
            const MBLevelMap<MAP, MEM>& map,
            MBIndex index) const;

        MBDisjointBoxLayout m_layout;
        std::set<Point> m_baseFootprint;
        Array<bool, DIM> m_copyAxes;
    };
    #include "implem/Proto_MBInterpLayoutImplem.H"

} //end namespace Proto
#endif // end include guard
