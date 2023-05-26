#pragma once
#ifndef _PROTO_MB_LEVEL_MAP_
#define _PROTO_MB_LEVEL_MAP_

#include "Proto_MBDisjointBoxLayout.H"
#include "Proto_MBLevelBoxData.H"
#include "Proto_Operator.H"

namespace Proto
{

template<MemType MEM = MEMTYPE_DEFAULT>
class MBLevelMap
{
    public:

    inline MBLevelMap(){}
    inline MBLevelMap(
            MBDisjointBoxLayout& a_layout,
            Array<Point, DIM+1>& a_ghost,
            int a_order = 4);

    inline void define(
            MBDisjointBoxLayout& a_layout,
            Array<Point, DIM+1>& a_ghost,
            int a_order = 4);

    /// Node Centered Map
    inline MBLevelBoxData<double, DIM, MEM, PR_NODE>& map();
    inline const MBLevelBoxData<double, DIM, MEM, PR_NODE>& map() const;
    
    /// Jacobian Access
    inline MBLevelBoxData<double, 1, MEM, PR_CELL>& jacobian();
    inline const MBLevelBoxData<double, 1, MEM, PR_CELL>& jacobian() const;

    /// User Defined Initialization
    inline virtual void init() {}

    /// Compute Map
    /**
     *  \param a_X:     Coordinates at nodes [output]
     *  \param a_J:     Cell averaged jacobian [output]
     *  \param a_index: MBIndex corresponding to the patch where the operator is being applied
     */
    inline virtual void apply(
            BoxData<double, DIM, MEM>& a_X,
            BoxData<double, 1, MEM>& a_J,
            const MBIndex& a_index);
    
    /// Compute Map (With Metrics)
    /**
     *  \param a_X:     Coordinates at nodes [output]
     *  \param a_J:     Cell averaged jacobian [output]
     *  \param a_NT:    Face averaged metric terms in each coordinate direction [output]
     *  \param a_index: MBIndex corresponding to the patch where the operator is being applied
     */
    inline virtual void apply(
            BoxData<double, DIM, MEM>& a_X,
            BoxData<double, 1, MEM>& a_J,
            FluxBoxData<double, DIM, MEM>& a_NT,
            const MBIndex& a_index);

    /// Get Layout
    inline const MBDisjointBoxLayout& layout() const;
    /// Get Mapped Grid Spacing
    inline const Array<double, DIM>& dx(int a_blockID) const;
    /// Get Mapped Coordinate Values
    /** 
     *  /param 
     */
    inline BoxData<double, DIM, MEM> X(const Box& a_box, const Array<double, DIM>& a_dx) const;

    /// Compute Mapped Point From Array
    /** a_X is an array in mapped space in the specified block
     *  The returned array is in physical space
     *
     */
    inline Array<double, DIM> operator()(
            const Array<double, DIM> a_X,
            unsigned int a_block) const;
    
    /// Compute Mapped Point From MBDataPoint
    inline Array<double, DIM> operator()(
            const MBDataPoint& a_pt,
            Point a_bound = Point::Zeros()) const;
    
    /// Compute Mapped Point From Point
    inline Array<double, DIM> operator()(
            const Point& a_cell,
            unsigned int a_block,
            Centering a_ctr) const;

    /// Compute Mapped Region
    inline BoxData<double, DIM> operator()(
            const Box& a_domain,
            unsigned int a_computeBlock,
            unsigned int a_domainBlock,
            Centering a_ctr = PR_NODE) const;
    
    /// Compute Mapped Region (No Reference Block)
    inline BoxData<double, DIM> operator()(
            const Box& a_domain,
            unsigned int a_block,
            Centering a_ctr = PR_NODE) const;
    
    // Find Cell
    inline bool contains(
            const MBIndex& a_patch, const Point& a_cell, const Array<double, DIM>& a_X) const;
    
    private:
        int                                         m_order;
        Array<Point, DIM+1>                         m_ghost;
        std::vector<Array<double, DIM>>             m_dx;
        MBLevelBoxData<double, DIM, MEM, PR_NODE>   m_X;
        MBLevelBoxData<double, 1, MEM, PR_CELL>     m_J;
        //MBLevelBoxData<double, DIM, MEM, PR_FACE>   m_NT;
        Stencil<double>                             m_c2c;
};

#include "implem/Proto_MBLevelMapImplem.H"
} // end namespace Proto
#endif // end include guard
