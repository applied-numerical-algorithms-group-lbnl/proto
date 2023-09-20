#pragma once
#ifndef _PROTO_MB_MAP_OP_
#define _PROTO_MB_MAP_OP_

#include "Proto_MBDisjointBoxLayout.H"
#include "Proto_MBLevelBoxData.H"
#include "Proto_Operator.H"

namespace Proto
{

class MBMapOp
{
    public:

    /// Trivial Constructor
    inline MBMapOp(){}

    /// Constructor
    inline MBMapOp(Array<double, DIM> a_dx, unsigned int a_block = 0);

    /// Lazy Constructor
    inline void define(Array<double, DIM> a_dx, unsigned int a_block = 0);

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
    template<MemType MEM>
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
    template<MemType MEM>
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
    template<MemType MEM>
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
    template<MemType MEM>
    inline BoxData<double, DIM, MEM> cellCentered(
            const Box& a_box,
            unsigned int a_computeBlock,
            unsigned int a_outBlock);
    
    /// Compute Cell Centered Coordinates
    /** Computes the cell-centered coodinate values at a specified MBDataPoint */
    inline Array<double, DIM> cellCentered(const MBDataPoint& a_point);
    
    /// Compute Cell Averaged Coordinates
    /** Computes the cell-averaged coodinate values at a specified MBDataPoint */
    inline Array<double, DIM> cellAveraged(const MBDataPoint& a_point);
    
    /// Compute Cell Centered Coordinates
    /** Returns the cell-averaged coordinates on a specified range.
     *
     *  \param a_box            Output range
     *  \param a_computeBlock   Block index associated with the coordinates being computed
     *  \param a_outBlock       Block index associated with the input range
     */
    template<MemType MEM>
    inline BoxData<double, DIM, MEM> cellAveraged(
            const Box& a_box,
            unsigned int a_computeBlock,
            unsigned int a_outBlock);
    
    /// Get Mapped Grid Spacing
    inline const Array<double, DIM>& dx(int a_blockID) const;

    /// Get Block
    inline unsigned int block() const;
    
    /// Get Mapped Coordinate Values
    /** Given a node-centered box and grid spacing, generate the node centered mapped
     *  coordinates. This assumes that the origin corresponds to the point (0,0,...,0) and
     *  that the input box already accounts for node centering (e.g. Box::grow(PR_NODE) is 
     *  assumed to already have been called.)
     *
     *  /param a_box     A box that already accounts for node centering
     *  /param a_dx      Grid spacing
     */
    template<MemType MEM>
    inline BoxData<double, DIM, MEM> X(const Box& a_box, const Array<double, DIM>& a_dx) const;
    
    private:
    
    std::vector<Array<double, DIM>> m_dx;   ///< Grid spacing in mapped space
    Stencil<double>                 m_c2c;  ///< Cached corners-to-cells Stencil
    unsigned int                    m_block;///< Block associated with this map (often unused)
};
#include "implem/Proto_MBMapOpImplem.H"
} //end Proto namespace
#endif //end include guard
