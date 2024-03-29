#pragma once
#ifndef _PROTO_MB_LEVEL_MAP_
#define _PROTO_MB_LEVEL_MAP_

#include "Proto_MBDisjointBoxLayout.H"
#include "Proto_MBLevelBoxData.H"
#include "Proto_Operator.H"

namespace Proto
{

/// Single Level Mapped Multiblock Map
/** This class provides an interface for constructing maps for use in single level
 *  mapped multiblock applications.
 *
 *  When a map is constructed, coordinate and Jacobian data are cached for all valid
 *  regions of the supplied MBDisjointBoxLayout as well as proscribed ghost regions.
 *  However, the map's apply function can be used to compute these quantities on any
 *  domain where the associated functions are defined. 
 *
 *  USAGE:
 *  To use this interface, create a class which
 *  inherits from MBLevelMap. The map itself is implemented by defining the one
 *  of the virtual "apply" functions (the user may supply either one or both
 *  depending on whether or not they need and/or would like to define how the metric
 *  tensor is computed). The user may also implement the virtual init function which will be
 *  called during construction of the Map. The init function is useful for caching constant
 *  data such as operators, stencils, or metrics in order to prevent computing these
 *  quantities each time apply is called.
 */
template<MemType MEM = MEMTYPE_DEFAULT>
class MBLevelMap
{
    public:

    /// Trivial Constructor
    inline MBLevelMap(){}
    
    /// Non-Trivial Constructor is explicitly removed
    /*
    inline MBLevelMap(
            MBDisjointBoxLayout& a_layout,
            Array<Point, DIM+1>& a_ghost);
    */
    
    /// Define
    /** Delayed construction. Define must be explicitly called by the user on an
     *  already declared map. This requirement ensures that the user defined
     *  initialization is called correctly. This restriction may be removed later.
     *
     *  \param a_layout     A mapped multiblock layout defining the domain
     *  \param a_ghost      An array of ghost data sizes. This ghost array is used
     *                      To define how much ghost data should be allocated for
     *                      cached arrays of the coordinates and Jacobian.
     *  \param a_block      The block id associated with this map. Often unused.
     */
    inline void define(
            const MBDisjointBoxLayout& a_layout,
            const Array<Point, DIM+1>& a_ghost,
            unsigned int         a_block = 0);

    /// Access Cached Coordinate Values
    /** Returns an MBLevelBoxData containing the cached coordinate data.
     *  The returned data has NODE centering.
     */
    inline MBLevelBoxData<double, DIM, MEM, PR_NODE>& map();
    
    /// Access Cached Coordinate Values (Const Overload)
    inline const MBLevelBoxData<double, DIM, MEM, PR_NODE>& map() const;
    
    /// Jacobian Access
    /** Returns an MBLevelBoxData containing the cached Jacobian data.
     *  The returned data has CELL centering and is assumed to be cell averaged.
     */
    inline MBLevelBoxData<double, 1, MEM, PR_CELL>& jacobian();
    
    /// Jacobian Access (Const Overload)
    inline const MBLevelBoxData<double, 1, MEM, PR_CELL>& jacobian() const;

    /// User Defined Initialization
    /** Users should override this function with any code they would like to run
     *  when the map is constructed.
     */
    inline virtual void init() {}

    /// Compute Map
    /** Users should override this function when defining their map if they do not want
     *  to explicitly define the metric tensor (NT). The domains of the coordinate and
     *  Jacobian quantities are assumed to be defined before they are input to this function.
     *  Note that the ghost values used to define the map during construction have absolutely
     *  no bearing on the valid domains for this function, so long as the relevant analytic
     *  functions are defined.
     *
     *  \param a_X:     Coordinates at nodes [input/output]
     *  \param a_J:     Cell averaged jacobian [input/output]
     *  \param a_index: MBIndex corresponding to the patch where the operator is being applied
     */
    inline virtual void apply(
            BoxData<double, DIM, MEM>& a_X,
            BoxData<double, 1, MEM>& a_J,
            unsigned int a_block);
    
    /// Compute Map (With Metrics)
    /** Users should override this function when defining their map if they want
     *  to explicitly define the metric tensor (NT). The domains of the coordinate and
     *  Jacobian quantities are assumed to be defined before they are input to this function. 
     *  Note that the ghost values used to define the map during construction have absolutely
     *  no bearing on the valid domains for this function, so long as the relevant analytic
     *  functions are defined.
     *
     *  \param a_X:     Coordinates at nodes [input/output]
     *  \param a_J:     Cell averaged jacobian [input/output]
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
     *  with a_index, computes X and J, and copies the result with rotation into the outputs.
     *
     *  Not recommended for public use.
     */
    inline void doApply(
            BoxData<double, DIM, MEM>& a_X,
            BoxData<double, 1, MEM>& a_J,
            unsigned int a_computeBlock,
            unsigned int a_outBlock);
   
    /// Compute Cell Centered Coordinates
    /** Returns the cell-centered coordinates on a specified range.
     *
     *  \param a_box            Output range
     *  \param a_computeBlock   Block index associated with the coordinates being computed
     *  \param a_outBlock       Block index associated with the input range
     */
    inline BoxData<double, DIM, MEM> cellCentered(
            const Box& a_box,
            unsigned int a_computeBlock,
            unsigned int a_outBlock);
    
    /// Compute Cell Centered Coordinates
    /** Returns the cell-averaged coordinates on a specified range.
     *
     *  \param a_box            Output range
     *  \param a_computeBlock   Block index associated with the coordinates being computed
     *  \param a_outBlock       Block index associated with the input range
     */
    inline BoxData<double, DIM, MEM> cellAveraged(
            const Box& a_box,
            unsigned int a_computeBlock,
            unsigned int a_outBlock);

    /// Compute Cell Centered Coordinates
    /** Computes the cell-centered coodinate values at a specified MBDataPoint */
    inline Array<double, DIM> cellCentered(const MBDataPoint& a_point);
    
    /// Compute Cell Averaged Coordinates
    /** Computes the cell-averaged coodinate values at a specified MBDataPoint */
    inline Array<double, DIM> cellAveraged(const MBDataPoint& a_point);

    /// Get Layout
    inline const MBDisjointBoxLayout& layout() const;
    
    /// Get Mapped Grid Spacing
    inline const Array<double, DIM>& dx(int a_blockID) const;

    /// Get Block
    inline unsigned int block() const { return m_block; }

    /// Get Mapped Coordinate Values
    /** Given a node-centered box and grid spacing, generate the node centered mapped
     *  coordinates. This assumes that the origin corresponds to the point (0,0,...,0) and
     *  that the input box already accounts for node centering (e.g. Box::grow(PR_NODE) is 
     *  assumed to already have been called.)
     *
     *  /param a_box     A box that already accounts for node centering
     *  /param a_dx      Grid spacing
     */
    inline BoxData<double, DIM, MEM> X(const Box& a_box, const Array<double, DIM>& a_dx) const;

    private:
        std::vector<Array<double, DIM>>             m_dx;   ///< Grid spacing in mapped space
        MBLevelBoxData<double, DIM, MEM, PR_NODE>   m_X;    ///< Cached coodinate values
        MBLevelBoxData<double, 1, MEM, PR_CELL>     m_J;    ///< Cached Jacobian values
        Stencil<double>                             m_c2c;  ///< Cached corners-to-cells Stencil
        unsigned int                                m_block;///< Block associated with this map (often unused)
};

#include "implem/Proto_MBLevelMapImplem.H"
} // end namespace Proto
#endif // end include guard
