#pragma once
#ifndef _PROTO_MB_LEVEL_MAP_
#define _PROTO_MB_LEVEL_MAP_

#include <optional>

#include "Proto_MBDisjointBoxLayout.H"
#include "Proto_MBLevelBoxData.H"
#include "Proto_MBMapOp.H"
#include "Proto_Operator.H"
//#include "Proto_CubedSphereKernels.H"

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
template<typename MAP, MemType MEM = MEMTYPE_DEFAULT>
class MBLevelMap
{
    public:

    /// Trivial Constructor
    inline MBLevelMap(){}
    
    /// Non-Trivial Constructor
    inline MBLevelMap(
            const MBDisjointBoxLayout& a_layout,
            const Array<Point, DIM+1>& a_ghost,
            unsigned int         a_focalBlock = 0);
    
    /// Non-Trivial Constructor
    inline MBLevelMap(
            const MBDisjointBoxLayout& a_layout,
            Point                      a_ghost,
            unsigned int         a_focalBlock = 0);
    
    /// Define
    /** Delayed construction.
     *
     *  \param a_layout     A mapped multiblock layout defining the domain
     *  \param a_ghost      An array of ghost data sizes. This ghost array is used
     *                      To define how much ghost data should be allocated for
     *                      cached arrays of the coordinates and Jacobian.
     *  \param a_focalBlock The "focus" block id associated with this map. Often unused.
     */
    inline void define(
            const MBDisjointBoxLayout& a_layout,
            const Array<Point, DIM+1>& a_ghost,
            unsigned int               a_focalBlock = 0);
    
    inline void define(
            const MBDisjointBoxLayout& a_layout,
            Point                      a_ghost,
            unsigned int               a_focalBlock = 0);

    inline void initialize() const;

    /// Access Cached Coordinate Values
    /** Returns an MBLevelBoxData containing the cached coordinate data.
     *  The returned data has NODE centering.
     */
    inline MBLevelBoxData<double, DIM, MEM, PR_NODE>& map();

    inline std::shared_ptr<MBLevelBoxData<double, DIM, MEM, PR_NODE>> mapData() {return m_X; }

    /// Access Cached Coordinate Values (Const Overload)
    inline const MBLevelBoxData<double, DIM, MEM, PR_NODE>& map() const;
    
    /// Jacobian Access
    /** Returns an MBLevelBoxData containing the cached Jacobian data.
     *  The returned data has CELL centering and is assumed to be cell averaged.
     */
    inline MBLevelBoxData<double, 1, MEM, PR_CELL>& jacobian();
    
    inline std::shared_ptr<MBLevelBoxData<double, 1, MEM, PR_CELL>> jacobianData() {return m_J; }
    /// Jacobian Access (Const Overload)
    inline const MBLevelBoxData<double, 1, MEM, PR_CELL>& jacobian() const;

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
    inline void apply(
            BoxData<double, DIM, MEM>& a_X,
            BoxData<double, 1, MEM>& a_J,
            unsigned int a_block) const;
    
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
    inline void apply(
            BoxData<double, DIM, MEM>& a_X,
            BoxData<double, 1, MEM>& a_J,
            FluxBoxData<double, DIM, MEM>& a_NT,
            unsigned int a_block) const;

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
            unsigned int a_outBlock) const;
   
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
            unsigned int a_outBlock) const;
   

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
            unsigned int a_outBlock) const;
    
    inline Array<double, DIM> mappedCoords(
        const Point& point,
        const BlockIndex& block) const;

    inline Array<double, DIM> mappedCoords(
        const Point& point,
        const BlockIndex& block,
        const Array<double, DIM>& offset) const;
    /// Compute Mapped Coordinates
    /** Returns coordinates in mapped space of a given data point. The Coordinates
     *  are of the cell center of the point plus the offset which is in units of the grid spacing. 
     */
//     inline Array<double, DIM> mappedCoords(
//             const MBDataPoint&          a_point,
//             const Array<double, DIM>&   a_offset) const;
    
//     /// Compute Mapped Coordinates (No Offset)
//     inline Array<double, DIM> mappedCoords(const MBDataPoint& a_point) const;

    /// Analytic Map (DataPoint Cell Center)
    inline Array<double, DIM> apply(const MBDataPoint& a_point) const;
    
    /// Analytic Map (DataPoint With Offset)
    inline Array<double, DIM> apply(
            const MBDataPoint&          a_point,
            const Array<double, DIM>&   a_offset) const;

    /// Compute Cell Centered Coordinates
    /** Computes the cell-centered coodinate values at a specified MBDataPoint */
    inline Array<double, DIM> cellCentered(const MBDataPoint& a_point) const;

    /// Compute Cell Centered Coordinates, analytic case.
    /** Computes the cell-centered coodinate valuesas a specified analytic function, */
    inline Array<double, DIM> analyticCellCentered(const MBDataPoint& a_point) const;
  
    /// Compute Cell Averaged Coordinates
    /** Computes the cell-averaged coodinate values at a specified MBDataPoint */
    inline Array<double, DIM> cellAveraged(const MBDataPoint& a_point) const;

    /// Get Layout
    inline const MBDisjointBoxLayout& layout() const;
    
    /// Get Mapped Grid Spacing
    inline const Array<double, DIM> dx(BlockIndex a_block) const;

    /// Get Operator
    inline MAP& op(BlockIndex a_block) const;

    /// Get Ghost Sizes
    inline Array<Point, DIM+1> ghost() const { return m_ghost; }

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

    inline MAP& operator[](BlockIndex block); 

    /// @brief Convert a Point from one block's index space to another by using forward / inverse mapping
    /// @param srcPoint Point in the source block's index space
    /// @param srcBlock Source block
    /// @param dstBlock Destination block
    /// @return 
    inline Point convertPoint(Point srcPoint, BlockIndex srcBlock, BlockIndex dstBlock) const;

    inline std::optional<std::pair<MBIndex, Point>> whichNeighborContains(const MBIndex& localIndex, const Point& localPoint, bool checkLocalInteriorGhostCells = true) const;

    private:
        Array<Point, DIM+1>             m_ghost;
        unsigned int                    m_focalBlock;
        Stencil<double>                 m_c2c;
        mutable std::vector<std::shared_ptr<MAP>>   m_ops;  ///< Map operators
        std::vector<Array<double, DIM>> m_dx;   ///< Mapped space grid spacing
        MBDisjointBoxLayout             m_layout;
        mutable bool                            m_initialized; 
        mutable std::shared_ptr<MBLevelBoxData<double, DIM, MEM, PR_NODE>>  m_X;    ///< Cached coodinate values
        mutable std::shared_ptr<MBLevelBoxData<double, 1, MEM, PR_CELL>>    m_J;    ///< Cached Jacobian values
        //mutable std::vector<std::unordered_map<MBPoint, MBPoint>>           m_pointMap; ///< Cached Point conversions
};

#include "implem/Proto_MBLevelMapImplem.H"
} // end namespace Proto
#endif // end include guard
