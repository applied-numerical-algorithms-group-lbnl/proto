#ifndef _PROTO_INTERP_STENCIL_1D_H_
#define _PROTO_INTERP_STENCIL_1D_H_
#include "Proto_Stencil.H"

namespace Proto {

/// Interpolation Stencil
/**
    \ingroup stencil_operations
    InterpStencil1D is a convenient abstraction which encapsulates the functionality of an operation in which
    the destination refinement is non-trivial. This kind of operation can be done with normal Stencil objects,
    however for a destination refinement of R, an array of R^DIM Stencils will be needed in general.
*/
template <class T>
class InterpStencil1D
{
public:
    /// Default Constructor
    inline InterpStencil1D() {}
    inline InterpStencil1D(unsigned int a_dir, unsigned int a_ratio) { }
    inline void define(unsigned int a_dir, unsigned int a_ratio);
    
    inline const Stencil<T>& get(unsigned int a_shift) const;
    inline void set(Stencil<T> a_S, unsigned int a_shift);

    /// Span
    /** \ingroup stencil_operations
        Computes the composite span of all Stencils in the InterpStencil1D
    */
    inline Box span() const;

    /// Ghost Size
    /**
        \ingroup stencil_operations
        Conservatively computes the amount of ghost data needed to apply this InterpStencil1D.
        May overestimate in the case of assymmetric InterpStencil1D footprints.
    */
    inline Point ghost() const;

    /// Lazy Apply
    /**
        \ingroup stencil_operations
        Gathers inputs for the composite Stencil operation and stores them in an intermediate
        LazyStencil to be consumed by a left hand side operator ("=", "|=", or "+=").
        The syntax is identical to that used by Stencil.

        \param a_src    Source data
        \param a_scale  (Optional) value to scale the output of *this by.
    */
    template<unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
    inline LazyStencil<T,C,MEMTYPE,D,E> operator()(
            const BoxData<T,C,MEMTYPE,D,E>& a_src,
            T                               a_scale = 1) const;
    
    /// Lazy Apply
    /**
        \ingroup stencil_operations
        Gathers inputs for the composite Stencil operation and stores them in an intermediate
        LazyStencil to be consumed by a left hand side operator ("=", "|=", or "+=").
        The syntax is identical to that used by Stencil.
        Explicitly defining a Box input is not recommended unless absolutely necessary.

        \param a_src    Source data
        \param a_box    Iteration box
        \param a_scale  (Optional) value to scale the output of *this by.
    */
    template<unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
    inline LazyStencil<T,C,MEMTYPE,D,E> operator()(
            const BoxData<T,C,MEMTYPE,D,E>& a_src,
            Box                             a_box,
            T                               a_scale = 1) const;

    /// Return Destination Refinement Ratio
    /**
        \ingroup stencil_operations
    */
    inline unsigned int ratio() const {return m_r;}
    /// Iterate Over Stencil Components
    /**
        \ingroup stencil_operations
    */
    inline typename std::vector<Stencil<T>>::iterator begin(){return m_s.begin();}
    /// Iterate Over Stencil Components
    /**
        \ingroup stencil_operations
    */
    inline typename std::vector<Stencil<T>>::iterator end(){return m_s.end();}

    inline void print() const;

#ifdef PR_OPS
//=======================================================================================
// InterpStencil1D Library ||
//=======================++
    #if 0
    /// Piecewise Constant Interpolation
    /**
        \ingroup stencil_library
        \param a_ratio  Anisotropic destination refinement
    */
    static inline InterpStencil1D<T> PiecewiseConstant(Point a_ratio);
    /// Piecewise Linear Interpolation
    /**
        \ingroup stencil_library
        \param a_ratio  Anisotropic destination refinement
    */
    static inline InterpStencil1D<T> PiecewiseLinear(Point a_ratio);
    /// Piecewise Linear Interpolation
    /**
        Computed based on tensor products of 1 dimensional interpolations.
        \ingroup stencil_library
        \param a_ratio  Isotropic destination refinement
    */
    static inline InterpStencil1D<T> Quadratic(int a_ratio);

    /// Build High-Order Interpolation Stencil (Simple)
    /**
        Builds a symmetric, <code>a_order</code> order accurate interpolation operator
        with an isotropic refinement ratio.
        
        \param a_order      Order of accuracy
        \param a_refRatio   Isotropic refinement ratio
    */
    static inline InterpStencil1D<T> Build(int a_order, int a_refRatio);
    
    /// Build High-Order Interpolation Stencil (Simple, Anisotropic)
    /**
        Builds a symmetric, <code>a_order</code> order accurate interpolation operator
        with an isotropic refinement ratio.
        
        \param a_order      Order of accuracy
        \param a_refRatio   Isotropic refinement ratio
    */
    static inline InterpStencil1D<T> Build(int a_order, Point a_refRatio);

    /// Build General High-Order Interpolation Stencil (General)
    /**
        Builds an InterpStencil1D with a custom footprint using a Least-Squares solve.

        \ingroup stencil_library
        \param a_shiftMax       Maximum Manhattan distance from the Stencil center
        \param a_shiftKernel    Maximum span of the InterpStencil1D
        \param a_order          Order of polynomial fit (e.g. 2 for quadratic)
        \param a_refRatio       Desired ratio of interpolation. Must be a power of 2
    */
    static inline InterpStencil1D<T> Build(
            int a_shiftMax,
            Box a_shiftKernel,
            int a_order,
            int a_refRatio);

    /// Build General High-Order Interpolation Stencil (General, Anisotropic)
    /**
        Builds an InterpStencil1D with a custom footprint using a Least-Squares solve.
        
        \ingroup stencil_library
        \param a_shiftMax       Maximum Manhattan distance from the Stencil center
        \param a_shiftKernel    Maximum span of the InterpStencil1D
        \param a_order          Order of polynomial fit (e.g. 2 for quadratic)
        \param a_refRatio       Desired ratio of interpolation. Must be a power of 2
    */
    static inline InterpStencil1D<T> Build(
            int a_shiftMax,
            Box a_shiftKernel,
            int a_order,
            Point a_refRatio);

    /// Build General High-Order Interpolation Stencil
    /**
        Builds an InterpStencil1D with a custom footprint using a Least-Squares solve.
        
        \ingroup stencil_library
        \param a_shifts         Desired shifts of the stencil.
        \param a_maxOrder       Highest desired order of Polynomial interpolant.
        \param a_refRatio       Desired ratio of interpolation. Must be a power of 2.
    */
    static inline InterpStencil1D<T> Build(
            std::vector<Point>& a_shifts,
            int                 a_maxOrder,
            Point               a_refRatio);
    
    /// Build General High-Order Interpolation Stencil
    /**
        Builds an InterpStencil1D with a custom footprint using a Least-Squares solve.
        Orders are specified as Point such that (1,2,3) is equivalent to x*y^2*z^3.

        \ingroup stencil_library
        \param a_shifts         Desired shifts of the stencil.
        \param a_orders         Polynomial exponents in vector form.
        \param a_refRatio       Desired ratio of interpolation. Must be a power of 2
    */
    static inline InterpStencil1D<T> Build(
            std::vector<Point>&       a_shifts,
            const std::vector<Point>& a_orders,
            Point                     a_refRatio);
#endif // end ifdef PR_OPS
#endif
private:
    unsigned int m_r;               ///< Destination refinement ratio
    unsigned int m_dir;              ///< Coordinate
    std::vector<Stencil<T>> m_s;    ///< Array of Stencils
};

#include "implem/Proto_InterpStencil1DImplem.H"
} // end namespace Proto
#endif // End Include Guard
