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
    
    /// Non-Trivial Constructor is explicitly removed
    /*
    inline MBLevelMap(
            MBDisjointBoxLayout& a_layout,
            Array<Point, DIM+1>& a_ghost,
            int a_order = 4);
    */
    
    inline void define(
            MBDisjointBoxLayout& a_layout,
            Array<Point, DIM+1>& a_ghost,
            unsigned int a_block = 0,
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
            unsigned int a_block);
    
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
            unsigned int a_block);

    /// Compute Map With Rotation
    /**
     *  Used internally for the case where X and J are defined on domains in 
     *  a different coordinate system from the block of the map computation (defined by a_index).
     *  This function creates a temporary version of X and J with domains in the block associated
     *  with a_index, computes X and J, and then copies the result with rotation into the outputs.
     */
    inline void doApply(
            BoxData<double, DIM, MEM>& a_X,
            BoxData<double, 1, MEM>& a_J,
            unsigned int a_computeBlock,
            unsigned int a_outBlock);
   
    /// Compute Cell Centered Coordinates
    inline BoxData<double, DIM, MEM> cellCentered(
            const Box& a_box,
            unsigned int a_computeBlock,
            unsigned int a_outBlock);
    
    /// Compute Cell Averaged Coordinates
    inline BoxData<double, DIM, MEM> cellAveraged(
            const Box& a_box,
            unsigned int a_computeBlock,
            unsigned int a_outBlock);

    inline Array<double, DIM> cellCentered(const MBDataPoint& a_point);
    
    inline Array<double, DIM> cellAveraged(const MBDataPoint& a_point);

    /// Get Layout
    inline const MBDisjointBoxLayout& layout() const;
    
    /// Get Mapped Grid Spacing
    inline const Array<double, DIM>& dx(int a_blockID) const;
    
    /// Get Mapped Coordinate Values
    /** Given a node-centered box and a grid spacing, generate the node centered mapped
     *  coordinates. This assumes that the origin corresponds to the point (0,0,...,0) and
     *  that the input box already accounts for node centering (e.g. Box::grow(PR_NODE) is 
     *  assumed to already have been called.)
     *
     *  /param a_box     A box that already accounts for node centering
     *  /param a_dx      Grid spacing
     */
    inline BoxData<double, DIM, MEM> X(const Box& a_box, const Array<double, DIM>& a_dx) const;

    /// Get Block
    inline unsigned int block() const { return m_block; }

    private:
        int                                         m_order;
        unsigned int                                m_block;
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
