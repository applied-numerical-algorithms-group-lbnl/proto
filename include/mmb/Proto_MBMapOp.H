#pragma once
#ifndef _PROTO_MB_MAP_OP_
#define _PROTO_MB_MAP_OP_

#include "Proto_MBDisjointBoxLayout.H"
#include "Proto_MBLevelBoxData.H"
#include "Proto_Operator.H"

namespace Proto
{

template<MemType MEM>
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
    
    /// Get Mapped Grid Spacing
    inline const Array<double, DIM>& dx() const;

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
    inline BoxData<double, DIM, MEM> X(const Box& a_box, const Array<double, DIM>& a_dx) const;
    
    private:
    
    Array<double, DIM> m_dx;   ///< Grid spacing in mapped space
    Stencil<double>    m_c2c;  ///< Cached corners-to-cells Stencil
    unsigned int       m_block;///< Block associated with this map (often unused)
};
#include "implem/Proto_MBMapOpImplem.H"
} //end Proto namespace
#endif //end include guard
