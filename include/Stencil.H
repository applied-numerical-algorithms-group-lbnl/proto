#ifndef _PROTO_STENCIL_H_
#define _PROTO_STENCIL_H_
#include "BoxData.H"

#include <vector>
#include <tuple>
#include <iostream>
#include <iomanip> //for pretty printing
namespace Proto {

//=======================================================================================
// SHIFT ||
//=======++ 

template <typename T> class Stencil; //forward declaration

/** @defgroup stencil_operations Stencil Operations*/
/*@{*/

/// Stencil Shift
/**
    \ingroup stencil_operations
    A shift is an alias for a Point which is used solely to provide
    a fluent syntax for creating Stencils. Refer to the documentation
    for Stencil for an illustrative example.
*/
class Shift {
public:
    
    ///////////////////////////////////////////////////////////////////////////////////////////////  
    /** @name Constructors */
    ///@{
  
    /// Default Constructor  
    Shift() : m_shift(Point::Zeros()){};

    /// Point Constructor
    /**
        \ingroup stencil_operations
        Builds a Shift from a Point
        
        \param a_pt   Point to base this Shift on
    */
    explicit Shift(const Point& a_pt){m_shift = a_pt;};

    /// Variadic Constructor
    /**
        \ingroup stencil_operations
        More or less identical to the variadic Point constructor.
        Builds a Shift using the first DIM arguments and ignores the rest.

        \param args   At least DIM int arguments, the first DIM of which define *this
    */
    template<typename... vals>
    inline explicit Shift(vals... args) : m_shift(args...) {}
    
    ///@}
    ///////////////////////////////////////////////////////////////////////////////////////////////  
    /** @name Methods */
    ///@{
    
    /// Basis Shift
    /**
        \ingroup stencil_operations
        Shortcut for Shift(Point::Basis(...))

        \param a_dir    Direction of basis vector in [0,DIM)
        \param a_scale  (OptionL) Amount to scale the basis vector by. (Default 1).
    */ 
    inline static Shift Basis(int a_dir, int a_scale=1)
    {return Shift(Point::Basis(a_dir, a_scale));};

    /// Zero Shift
    /**
        \ingroup stencil_operations
        Shortcut for Shift(Point::Zeros())
    */ 
    inline static Shift Zeros(){return Shift(Point::Zeros());};
    
    /// Unit Shift
    /**
        \ingroup stencil_operations
        Shortcut for Shift(Point::Ones(...))

        \param a_scale    (Optional) Used to scale the default output of (1,1,...,1).
    */ 
    inline static Shift Ones(int a_scale=1){return Shift(Point::Ones(a_scale));};

    /// Get Shift Point
    inline Point& shift(){return m_shift;}

    /// Scalar Multiplication
    /**
        \ingroup stencil_operations
        Generates a Stencil<T> from the product of a scalar T coefficient and a Shift.
        This operator is what allows for a fluid Stencil construction syntax:

        Example:
        @code
        // DIM == 1
        Stencil<double> = 1.0*Shift(-1) - 2.0*Shift(0) + 1.0*Shift(1);
        @endcode
    */
    template<typename T>
    inline Stencil<T> operator*(T a_coef) const
    {
        return Stencil<T>(*this,a_coef);
    }
    
    /// Convolution
    /**
        \ingroup stencil_operations
        The product of two Shift objects is defined as their sum. 
    */
    inline Shift operator*(const Shift& a_shift) const
    {
        return Shift(m_shift + a_shift.m_shift);
    }

    /// Componentwise Access
    inline int& operator[](int a_dir)
    {
      PROTO_ASSERT(((a_dir >= 0) && (a_dir < DIM)),
        "Shift::operator[](int a_dir) invalid for a_dir = %i. a_dir must be in [0,DIM=%i)",
        a_dir, DIM);
        return m_shift[a_dir];
    }
    */
    ///@}
private:
    Point m_shift;
};

//=======================================================================================
// LAZYSTENCIL ||
//=============++ 

// forward declaration
template <typename T, unsigned int C, unsigned char D, unsigned char E>
class BoxData;

/// An Unevaluated Stencil Operation
/**
    \ingroup stencil_operations
    LazyStencil is an intermediate structure that holds the intermediate data for a 
    Stencil operation. 

    LazyStencil is not explicitly part of the user interface, and is only public by 
    virtue of necessity. 
*/
template <typename T, unsigned int C, unsigned char D, unsigned char E>
struct LazyStencil {
    
    inline LazyStencil() {};
    inline LazyStencil(Stencil<T> a_stencil,
                       const BoxData<T,C,D,E>* a_src,
                       Bx a_box, T a_scale);
    
    inline void apply(BoxData<T,C,D,E>& a_dest, bool a_overwrite);

    inline int size(){return m_stencil.size();}
    
    Bx m_range;

    std::vector<Stencil<T>> m_stencil;
    std::vector<BoxData<T, C, D, E>*> m_src;
    std::vector<Bx> m_box;
    std::vector<T> m_scale;
};

//=======================================================================================
// STENCIL ||
//=========++ 

/// A Linear Stencil Operation
/**
    \ingroup stencil_operations
    Encapsulates a linear stencil operation where coefficients are of type T.
    Stencil objects are built and used in a way that conforms to their nature as operators.
    For illustrative usage examples, refer to the following code snippets:

    Examples:
    Build a Stencil from Shifts and coefficients:
    \snippet Snippets.cpp proto_stencil_build
    Apply a Stencil with no Source / Dest Refinement to a BoxData:
    \snippet Snippets.cpp proto_stencil_apply

    The above examples illustrate the usage of Stencil to do computations in which the source and
    destination arrays are at the same refinement. Stencil is also capable of "prolong" and "restrict"
    type operations in which the destination and source respectively are refined. This functionality is 
    useful for operations such as averaging and interpolation, or for algorithms like Red-Black Gauss Seidel
    Iteration in which it is necessary to evaluate a Stencil on "every other cell". 

    To facilitate these more exotic operations, the Stencil API allows the user to designate a source and/or destination
    refinement ratio. If these values are different from (1,...,1), then input Bx object will be interpreted as an
    "index range" instead of a physical Bx domain. The following code snippets illustrate some examples of this functionality.

    Examples:
    Non-Trivial Source Refinement Ratio
    \snippet Snippets.cpp proto_stencil_average
    Non-Trivial Destination Refinement Ratio
    \snippet Snippets.cpp proto_stencil_dest_refine

    In the case of non-trivial destination refinement, an array of refRatio^DIM Stencils is often needed to fully populate the output.
    Proto provides a convenience structure called InterpStencil which is designed to mitigate this form of pedantry. See the associated
    documentation for additional information and example usage.
*/
template <typename T>
class Stencil {

public:

    template <typename TT, unsigned int C, unsigned char D, unsigned char E>
    friend class LazyStencil;

    template <typename TT, unsigned int C, unsigned char D, unsigned char E>
    friend class BoxData;

    template <typename TT, unsigned int C, unsigned char D, unsigned char E>
    friend BoxData<TT,C,D,E>& operator|=(BoxData<TT,C,D,E>& a_dest,
                                         LazyStencil<TT,C,D,E>&& a_op);

    template <typename TT, unsigned int C, unsigned char D, unsigned char E>
    friend BoxData<TT,C,D,E>& operator+=(BoxData<TT,C,D,E>& a_dest,
                                         LazyStencil<TT,C,D,E>&& a_op);
    
    ///////////////////////////////////////////////////////////////////////////////////////////////  
    /** @name Constructors */
    ///@{

    /// Default Constructor
    Stencil();
    
    /// General Constructor
    /**
        \ingroup stencil_operations
        Creates a Stencil with a single shift and coefficent.

        Not recommended for public use; see the Stencil class documentation for examples of how to 
        build a Stencil with Shift - coefficient syntax.

        \param a_shift        Shift of this operation
        \param a_coef         Coefficient of this operation
        \param a_destRefratio (Optional) Destination refinement ratio. Defaults to (1,...,1) (no refinement)
        \param a_destShift    (Optional) Destination shift. Defaults to (0,...,0) (no shift)
        \param a_srcRefratio  (Optional) Source refinement ratio. Defaults to (1,...,1) (no refinement)
    */
    Stencil(Shift   a_shift,
            T       a_coef,
            Point   a_destRefratio  = Point::Ones(),
            Point   a_destShift     = Point::Zeros(),
            Point   a_srcRefratio   = Point::Ones());

    ///@}

    ///////////////////////////////////////////////////////////////////////////////////////////////  
    /** @name Operators */
    ///@{

    /// Stencil Composition
    /**
        \ingroup stencil_operations
        The product of two Stencils is defined as the composition of the two.

        \param a_stencil  Another Stencil
    */
    Stencil<T> operator*(const Stencil<T>& a_stencil) const;
    
    /// Scalar Multiplication
    /**
        \ingroup stencil_operations
        The product of a Stencil S and a coefficient v of type T results in the scaling of all
        the coefficients of S by v.

        \param a_coef   Scaling coefficient
    */
    Stencil<T> operator*(const T a_coef) const;
 
    /// In Place Stencil Composition
    /**
        \ingroup stencil_operations
        \param a_stencil  Another Stencil
    */
    void operator*=(const Stencil<T>& a_stencil);

    /// In Place Scalar Multiplication
    void operator*=(const T a_coef);

    /// Stencil Addition
    /**
        \ingroup stencil_operations
        Adding two Stencils results in a new Stencil which is the union of the coefficent-Shift pairs of the inputs.
        If the two input Stencils share a common Shift, the associated coefficients will be added together.

        \param a_stencil  Another Stencil
    */
    Stencil<T> operator+(const Stencil<T>& a_stencil) const;
    
    /// Stencil Subtraction (Convenience)
    /**
        \ingroup stencil_operations
        Equivalent to adding <code> a_stencil*(-1) </code>

        \param a_stencil  Another Stencil
    */
    Stencil<T> operator-(const Stencil<T>& a_stencil) const;

    /// In Place Stencil Addition
    /**
        \ingroup stencil_operations
        \param a_stencil  Another Stencil
    */
    void operator+=(const Stencil<T>& a_stencil);
 
    /// In Place Stencil Subtraction
    /**
        \ingroup stencil_operations
        \param a_stencil  Another Stencil
    */
    void operator-=(const Stencil<T>& a_stencil);
 
    ///@}

    ///////////////////////////////////////////////////////////////////////////////////////////////  
    /** @name Accessors and Queries*/
    ///@{
    /// Equality Operator
    /**
        \ingroup stencil_operations
        Equality between Stencils is determined by value
        
        \param a_stencil  Another Stencil
    */
    inline bool operator==(Stencil<T>& a_stencil) const;
    
    /// Inquality Operator
    /**
        \ingroup stencil_operations
        \param a_stencil  Another Stencil
    */
    inline bool operator!=(Stencil<T>& a_stencil) const {return !(*this == a_stencil);}
    
    /// Get Vector of Coefficients
    /**
        \ingroup stencil_operations
        Returned vector is read-only. Ordering corresponds to the output of Stencil::offsets()
    */
    inline const std::vector<T>& coefs() const {return m_coefs;};

    /// Get Vector of Offsets
    /**
        \ingroup stencil_operations
        Returned vector is read-only. Ordering corresponds to the output of Stencil::coefs()
    */
    inline const std::vector<Point>& offsets() const {return m_offsets;};

    /// Size
    /**
        \ingroup stencil_operations
        Defined as the number of coefficient-offset pairs in the stencil
    */
    inline unsigned int size() const {return m_coefs.size();}

    /// Span
    /**
        \ingroup stencil_operations
        Returns a Box which bounds all offsets in *this. Useful for automated ghost-cell checks.

        Example:
        @code
        //DIM = 2
        using namespace Proto;
        Stencil<T> S = 1.0*Shift(0,-1) + 
                       2.0*Shift(0,2) + 
                       3.0*Shift(1,3);
        std::cout << S.span() << std::endl; //prints [(-1,0), (2,3)] 
        @endcode
    */
    inline Bx span() const {return m_span;};

    /// Ghost
    /**
        \ingroup stencil_operations
        Returns the necessary ghost region width in the form of a Point.
        If a Stencil requires a different number of ghost cells in the "upper" direction than the "lower" direction,
        the larger value will be used.
    */
    inline Point ghost() const
    {
        Point ghost;
        for (int ii = 0; ii < DIM; ii++)
        {
            ghost[ii] = std::max(-m_span.low()[ii], m_span.high()[ii]);
        }
        return ghost;
    }

    /// Get Source Refinement Ratio
    /**
        \ingroup stencil_operations
        returned reference is mutable. Use this function to set the source refinement ratio of custom Stencils
    */
    inline Point& srcRatio(){return m_srcRefratio;};
    
    /// Get Source Refinement Ratio (Const)
    inline const Point& srcRatio() const {return m_srcRefratio;};
    
    /// Get Destination Refinement Ratio
    /**
        \ingroup stencil_operations
        returned reference is mutable. Use this function to set the destination refinement ratio of custom Stencils
    */
    inline Point& destRatio(){return m_destRefratio;};
    
    /// Get Destination Refinement Ratio (Const)
    inline const Point& destRatio() const {return m_destRefratio;};
    
    /// Get Destination Shift
    /**
        returned reference is mutable. Use this function to set the destination shift of custom Stencils
    */
    inline Point& destShift(){return m_destShift;};
    
    /// Get Destination Shift (Const)
    /**
        \ingroup stencil_operations
    */
    inline const Point& destShift() const {return m_destShift;};
    
    ///@}

    ///////////////////////////////////////////////////////////////////////////////////////////////  
    /** @name Stencil Binding and Application */
    ///@{
    
    /// Operate on BoxData
    /**
        \ingroup stencil_operations
        Operate *this on a BoxData. This function works in tandem with the namespace-defined
        operators += and |=. See the example in the description of Stencil.

        \param a_scr    Source BoxData
        \param a_scale  (Optional) Scale the output of the Stencil operation. Useful for eliding calls to BoxData::operator*=(...) and similar.
    */
    template<unsigned int C, unsigned char D, unsigned char E>
    inline LazyStencil <T,C,D,E>
    operator()(const BoxData<T,C,D,E>&  a_src,
               T                        a_scale = 1) const 
    {
        //Bx b = (this->range(a_src.box())).taperCoarsen(m_destRefratio);
        Bx b = this->indexRange(a_src.box());
        return LazyStencil<T,C,D,E>(*this, &a_src, b, a_scale);
    }
    /// Operate on BoxData (Overload with Bx Input)
    /**
        \ingroup stencil_operations
        Operate *this on a BoxData. This function works in tandem with the namespace-defined
        operators += and |=. See the example in the description of Stencil.

        \param a_scr    Source BoxData
        \param a_box    Confinement Bx. Must be a subset of the allowable range of *this
        \param a_scale  (Optional) Scale the output of the Stencil operation. Useful for eliding calls to BoxData::operator*=(...) and similar.
    */
    template<unsigned int C, unsigned char D, unsigned char E>
    inline LazyStencil <T,C,D,E>
    operator()(const BoxData<T,C,D,E>&  a_src,
               Bx                       a_box,
               T                        a_scale = 1) const 
    {
    #ifndef NDEBUG
        Bx b = this->indexRange(a_src.box());
        PROTO_ASSERT(b.contains(a_box),
        "Stencil::operator() invalid. User provided a Bx larger than the valid range of the input data.");
    #endif
        return LazyStencil<T, C, D, E>(*this, &a_src, a_box, a_scale);
    }
    ///@}
    
    ///////////////////////////////////////////////////////////////////////////////////////////////  
    /** @name Utility */
    ///@{
    
    /// Transpose Stencil
    /**
        \ingroup stencil_operations
        Transposes *this across two directions in [0,DIM). After transpose, coefficients associated with
        the offset (a,...,b,...) will be associated instead with (b,...,a,...)
    */
    inline void transpose(unsigned char a, unsigned char b)
    {
        PROTO_ASSERT((a < DIM) && (b < DIM),
        "Stencil::transpose(int a, int b) invalid. Both a and b must be in [0,DIM = %i)",DIM);
        #if DIM < 2
        return;
        #else
        if (a == b) {return;}
        for (int l = 0; l < size(); l++)
        {
            std::swap(m_offsets[l][a], m_offsets[l][b]);
        }
        #endif
    }

    /// Get Max Index Range
    /**
        \ingroup stencil_operations
        Given a domain, compute the largest possible iteration Box, taking source/ destination refinement into account.
        For Stencils without source or destination refinement, this function is identical to Stencil::range.

        \param a_domain     A computational domain
    */
    inline Bx indexRange(Bx a_domain) const
    {
      Bx indexRange = Bx(a_domain.low() - m_span.low(),
                         a_domain.high() - m_span.high());
      indexRange = indexRange.taperCoarsen(m_srcRefratio);
      return indexRange;
    }
    
    /// Get Max Index Domain
    /**
        \ingroup stencil_operations
        Given a domain, compute the largest possible iteration Box, taking source / destination refinement into account.
        The output of this function is always a valid input for a Stencil operation.
        For Stencils without source or destination refinement, this function is identical to Stencil::domain.

        \param a_range      A computational range
    */
    inline Bx indexDomain(Bx a_range) const
    {
      Bx indexDomain = a_range.shift(-m_destShift);
      indexDomain = indexDomain.taperCoarsen(m_destRefratio);
      return indexDomain;
    }

    /// Get Max Range Box
    /**
        \ingroup stencil_operations
        Given a domain, compute the largest associated physical range, taking refinement AND destination shifting into account.
        The output of this function is NOT a valid input for a Stencil operation when refinement is present.
        This function is best used for defining output BoxData when the input domain is known.

        \param a_domain     A computational domain
    */
    inline Bx range(Bx a_domain) const
    {
        Bx range;
        range = indexRange(a_domain);
        range = Bx(range.low()*m_destRefratio, range.high()*m_destRefratio);
        range = range.shift(m_destShift);
        return range;
    }
    
    /// Get Min Domain Bx
    /**
        \ingroup stencil_operations
        Given a range, compute the smallest associated physical domain, taking refinement AND destination shifting into account.
        The output of this function is NOT a valid input for a Stencil operation when refinement is present.
        Thisfunction is best used for defining input BoxData when the output range is known.

        \param a_range     A computational range
    */
    inline Bx domain(Bx a_range) const
    {
        Bx domain = indexDomain(a_range);
        domain = Bx(domain.low()*m_srcRefratio, domain.high()*m_srcRefratio);
        domain = Bx(domain.low() + m_span.low(), domain.high()+m_span.high());
        return domain;
    }

    /// Print
    /**
        \ingroup stencil_operations
        Print *this to the command line. Useful for debugging.
    */
    inline void print() const
    {
        using ::std::cout;
        using ::std::endl;
        
        cout << "Coeficients and Shifts: { ";
        for (int ii = 0; ii < m_coefs.size(); ii++) {
            if (ii == 0){cout << endl;}
            cout << std::setfill(' ') << std::setw(8) << m_coefs[ii];
            cout << ": " << m_offsets[ii] << endl;
        }
        cout << "}" << endl;
        cout << "srcRefratio:\t" << m_srcRefratio << endl;
        cout << "destRefratio:\t" << m_destRefratio << endl;
        cout << "destShift:\t" << m_destShift << endl;
        cout << "span:\t\t" << m_span << endl;
        cout << endl;
    }
    ///@}

    // End of Stencil Operations Doxygen Module
    /*@}*/
    
    /** @defgroup stencil_library Stencil Library*/
    /*@{*/
    
    ///////////////////////////////////////////////////////////////////////////////////////////////  
    /** @name Stencil Library */
    ///@{
    
    /// Stencil Library: Derivative
    /**
        Built in implementation of compact differentiation stencils. 
        Includes derivatives of order n >= 1 and accuracy m >= 2 where n + m <= 14.

        \param a_n      Degree of derivative (e.g. first derivative, third deriative etc.)
        \param a_dir    (Optional) Direciton of derivative. An unsigned int in [0,DIM). Default 0
        \param a_order  (Optional) Order of accuracy. An unsigned int >= 2. Maximum accuracy depends on a_n. Default 2.
    */
    static Stencil<T> Derivative(int a_n, int a_dir = 0, int a_order = 2);
    
    /// Stencil Library: Laplacian
    /**
        Built in implementation of the 2nd order 2*DIM + 1 point Laplace operator.
    */
    static Stencil<T> Laplacian();
    #if DIM == 2
    /// 9 Point Laplacian
    static Stencil<T> Laplacian_9();
    #elif DIM == 3
    /// 27 Point Laplacian
    static Stencil<T> Laplacian_27();
    #endif
    /// Stencil Library: Perpendicular Laplacian
    /**
        Built in implementation of Laplacian perpendicular to direction dir
        
        \param a_dir    Normal direction
        \param a_order  (Optional) Order of accuracy (default: 2, supported: 2)
    */
    static Stencil<T> LaplacianFace(int a_dir, int a_order = 2);
    
    /// Stencil Library: Edge to Cell Convolution
    /**
        Interpolates values averaged at cell i to exact values at the cell boundary i+1/2.
        \param a_dir    Axis to interpolate along
        \param a_order  (Optional) Order of accuracy (default: 4, supported: 4)
    */
    static Stencil<T> EdgeToCell(int a_dir, int a_order = 4);

    /// Stencil Library: Downwind Cell to Edge Deconvolution 
    /**
        Interpolates values averaged at cell i to exact values at cell boundary i-1/2
        using downwind data.

        \param a_dir    Axis to interpolate along
        \param a_order  (Optional) Order of accuracy (default: 5, supported: 5)
    */
    static Stencil<T> CellToEdgeL(int a_dir, int a_order = 5);
    
    /// Stencil Library: Upwind Cell to Edge Deconvolution 
    /**
        Interpolates values averaged at cell i to exact values at cell boundary i-1/2
        using upwind data.

        \param a_dir    Axis to interpolate along
        \param a_order  Order of accuracy (default: 5, supported: 5)
    */
    static Stencil<T> CellToEdgeH(int a_dir, int a_order = 5);
    
    /// Stencil Library: Simple Average
    /**
        Averages data from a refined grid onto a coarsened grid.
        Refinement is assumed to be isotropic.

        \param a_refRatio   (Optional) Refinement of source data relative to destination (default: 2)
    */
    static Stencil<T> AvgDown(int a_refRatio = 2);
    
    /// Stencil Library: Flux Divergence
    /**
        Simple flux differencing stencil: OUT(i) = IN(i+1) - IN(i)
        Assumes the lower flux is stored in cell i.

        \param a_dir    Axis to difference along.
    */
    static Stencil<T> FluxDivergence(int a_dir);
    ///@}
    
    // End of Stencil Operations Doxygen Module
    /*@}*/
    
    /// Apply Stencil Helper function
    /**
        \ingroup stencil_operations
        Manually appy *this to a source and destination BoxData in a functional programming style.
        Not recommended for public use, but useful for debugging in some cases.

        \param a_src        Source data
        \param a_dest       Output data
        \param a_bx         Iteration box for this computation. Can be created with indexDomain or indexRange functions.
        \param a_replace    (Optional) If true, zero out the data in a_dest and replace with the computation result. If true, solution will be added to existing data. Default false.
        \param a_scale      (Optional) Scale the computation by some value. 
    */
    template<unsigned int C, unsigned char D, unsigned char E>   
    void apply(const BoxData<T,C,D,E>&   a_src,
               BoxData<T,C,D,E>&         a_dest,
               const Bx&                 a_bx,
               bool                      a_initToZero = false,
               const T                   a_scale = 1) const;


    template<unsigned int C, unsigned char D, unsigned char E>   
    void  getIntegerOffsets(std::vector<int>      &  a_srcOffsets,
                            std::vector<int>      &  a_dstOffsets,
                            int                   &  a_bxpts,
                            Bx              &  a_bx,
                            BoxData<T,C,D,E>&  a_src,
                            BoxData<T,C,D,E>&  a_dst) const;

#ifdef PROTO_CUDA

    template<unsigned int C, unsigned char D, unsigned char E>   
    void cudaApply( BoxData<T,C,D,E>&  a_src,
                    BoxData<T,C,D,E>&  a_dst,
                    Bx&                a_box,
                    bool               a_initToZero = false,
                    T                  a_scale = 1) ;                      

#endif
private:
    /// Augment Stencil
    /**
        Helper function that encapsulates all of the proper checks needed
        to add a coefficient-offset pair to a Stencil. Any new function 
        that needs to add data to an existing Stencil should call this.
    */
    void augment(T a_coef, Point a_offset);
    
    std::vector<T>      m_coefs;        ///< Coefficients of the Stencil.
    std::vector<Point>  m_offsets;      ///< Offsets associated with the Stencil.
    Point               m_srcRefratio;  ///< Refinement of source data
    Point               m_destRefratio; ///< Refinement of destination data
    Point               m_destShift;    ///< Output shift in (refined) destination data. Meaningless without destination refinement
    Bx                  m_span;         ///< Bounding Box defining the largest offsets of *this. Defines Stencil's spatial footprint
};

///////////////////////////////////////////////////////////////////////////////////////////////  
/** @name Non-Member Functions */
///@{

/// Coefficient Shift Product "Constructor"
/**
    \ingroup stencil_operations
    Syntactical sugar that allows a Stencil to be constructed with coef*Shift(Point) syntax
  
    Example:
    @code
    Stencil<double> S = 3.7*Shift(Point::Basis(0));
    @endcode
*/
template <typename T>
inline Stencil<T> operator*(T a_coef, Shift a_shift)
{
    return Stencil<T>(a_shift, a_coef);
}

/// Scalar Multiplication of Stencil Coefficients
/**
    \ingroup stencil_operations
    Allows for pre multiplication by a T scalar
*/
template <typename T>
inline Stencil<T> operator*(T a_coef, const Stencil<T> a_stencil)
{
    return a_stencil*a_coef;
}

/// Application by Replacement
/**
    \ingroup stencil_operations
    Applies a LazyStencil and replaces a subset of a_dest with the result:
    
    Usage:
    @code
    Stencil<double> S;
    BoxData<double> Src;
    Bx computeBox;
    // define S, computeBox and initialize Src
    BoxData<double> Dest |= S(Src, computeBox);
    @endcode
    
    See the main documentation for Stencil for additional examples
    
    \param a_dest   Destination array
    \param a_op     Uncomputed Stencil operation
*/
template <typename T, unsigned int C, unsigned char D, unsigned char E>
inline BoxData<T,C,D,E>& operator|=(BoxData<T,C,D,E>& a_dest, LazyStencil<T,C,D,E>&& a_op)
{
    if (a_op.size() > 0)
    {
        for (auto iter = a_op.m_src.begin(); iter != a_op.m_src.end(); ++iter)
        {
            PROTO_ASSERT(!(*iter)->isAlias(a_dest),
            "BoxData::operator|=(BoxData) invalid. Stencil source and destination data are aliased");
        }
        a_op.apply(a_dest, true);
    }
    return a_dest; 
}
 
/// Application by Increment
/**
    \ingroup stencil_operations
    Applies a LazyStencil and replaces a subset of a_dest with the result:
    
    Usage:
    @code
    Stencil<double> S;
    BoxData<double> Src;
    Bx computeBox;
    // define S, computeBox and initialize Src
    BoxData<double> Dest += S(Src, computeBox);
    @endcode
    
    See the main documentation for Stencil for additional examples
    
    \param a_dest   Destination array
    \param a_op     Uncomputed Stencil operation
*/
template <class T, unsigned int C, unsigned char D, unsigned char E>
BoxData<T,C,D,E>& operator+=(BoxData<T,C,D,E>& a_dest, LazyStencil<T,C,D,E>&& a_op)
{
    if (a_op.size() > 0)
    {
        for (auto iter = a_op.m_src.begin(); iter != a_op.m_src.end(); ++iter)
        {
            PROTO_ASSERT(!(*iter)->isAlias(a_dest),
            "BoxData::operator|=(BoxData) invalid. Stencil source and destination data are aliased");
        }
        a_op.apply(a_dest,false); 
    }
    return a_dest;
}

///@}

//=======================================================================================
// InterpStencil ||
//===============++ 

/// Interpolation Stencil
/**
    \ingroup stencil_operations
    InterpStencil is a convenient abstraction which encapsulates the functionality of an operation in which
    the destination refinement is non-trivial. This kind of operation can be done with normal Stencil objects,
    however for a destination refinement of R, an array of R^DIM Stencils will be needed in general. 
*/

template <class T>
class InterpStencil
{
public:
    /// Default Constructor
    inline InterpStencil()
    {
        m_closed = false;
    }
    
    /// Isotropic Constructor
    /**
        \ingroup stencil_operations
        Builds an InterpStencil assuming a destination refinement of <code> a_ratio </code> in all directions.

        \param a_ratio  destination refinement
    */
    inline InterpStencil(int a_ratio)
    {
        define(Point::Ones(a_ratio));
    }
    
    /// Anisotropic Constructor
    /**
        \ingroup stencil_operations
        Builds an InterpStencil with a destination refinement of a_ratio;

        \param a_ratio  destination refinement
    */
    inline InterpStencil(Point a_ratio)
    {
        define(a_ratio);
    }
    
    
    /// Define
    /**
        \ingroup stencil_operations
        Builds an InterpStencil with a destination refinement a_ratio.

        \param a_ratio  destination refinement
    */
    inline void define(Point a_ratio)
    {
        for (int ii = 0; ii < DIM; ii++)
        {
            PROTO_ASSERT(a_ratio[ii] > 0,
            "InterpStencil(Point ratio) invalid. All ratios must be 1 or greater");
        }
        m_k = Bx(a_ratio - Point::Ones());
        m_s.resize(m_k.size());
        m_closed = false;
    }

    inline void close()
    {
        for (auto iter = m_k.begin(); iter != m_k.end(); ++iter)
        {
            Stencil<T>& S = (*this)(*iter);
            S.destRatio() = m_k.high() + Point::Ones();
            S.destShift() = *iter;
        }
        m_closed = true;
    }

    /// Get Read-Write Stencil Subcomponent
    /**
        \ingroup stencil_operations
        Returns a reference to the Stencil associated with the destShift a_p. 

        \param a_p  A destination shift in [(0,...,0), (destRatio-1,...,destRatio-1)]
    */
    inline Stencil<T>& operator()(Point a_p)
    {
        PROTO_ASSERT(!m_closed,
        "Components of InterpStencil are read-only once closed. Use InterpStencil::get instead.");
        PROTO_ASSERT(m_k.contains(a_p),
        "InterpStencil::operator()(Point p) invalid.\
         p is not a member of the InterStencil kernel");
        
        return m_s[m_k.index(a_p)];
    }
    
    /// Get Read-Only Stencil Subcomponent
    /**
        \ingroup stencil_operations
        Returns a const reference to the Stencil associated with the destShift a_p. 

        \param a_p  A destination shift in [(0,...,0), (destRatio-1,...,destRatio-1)]
    */
    inline const Stencil<T>& get(Point a_p)
    {
        PROTO_ASSERT(m_k.contains(a_p),
        "InterpStencil::operator()(Point p) invalid.\
         p is not a member of the InterStencil kernel");

        return m_s[m_k.index(a_p)];
    }

    /// Lazy Apply
    /**
        \ingroup stencil_operations
        Gathers inputs for the composite Stencil operation and stores them in an intermediate LazyStencil
        to be consumed by a left hand side operator ("=", "|=", or "+="). The syntax is intentionally identical to that used by Stencil.
        

        \param a_src    Source data
        \param a_scale  (Optional) value to scale the output of *this by. 
    */
    template<unsigned int C, unsigned char D, unsigned char E>
    inline LazyStencil<T,C,D,E> operator()(const BoxData<T,C,D,E>& a_src,
                                           T                       a_scale = 1) const;
    /// Lazy Apply
    /**
        \ingroup stencil_operations
        Gathers inputs for the composite Stencil operation and stores them in an intermediate LazyStencil
        to be consumed by a left hand side operator ("=", "|=", or "+="). The syntax is intentionally identical to that used by Stencil.
        Explicitly defining a Bx input is not recommended unless absolutely necessary.

        \param a_src    Source data
        \param a_box    Iteration box
        \param a_scale  (Optional) value to scale the output of *this by. 
    */
    template<unsigned int C, unsigned char D, unsigned char E>
    inline LazyStencil<T,C,D,E> operator()(const BoxData<T,C,D,E>& a_src,
                                           Bx                      a_box,
                                           T                       a_scale = 1) const;
    /// Piecewise Constant Interpolation
    /**
        \ingroup stencil_library
        \param a_ratio  Anisotropic destination refinement
    */    
    static inline InterpStencil<T> PiecewiseConstant(Point a_ratio);
    /// Piecewise Linear Interpolation
    /**
        \ingroup stencil_library
        \param a_ratio  Anisotropic destination refinement
    */    
    static inline InterpStencil<T> PiecewiseLinear(Point a_ratio);
    /// Piecewise Linear Interpolation
    /**
        Computed based on tensor products of 1 dimensional interpolations.
        \ingroup stencil_library
        \param a_ratio  Isotropic destination refinement
    */    
    static inline InterpStencil<T> Quadratic(int a_ratio);
   
    /// Return Shift Kernel 
    /** 
        \ingroup stencil_operations
    */
    inline Bx kernel() const {return m_k;}
    /// Return Destination Refinement Ratio 
    /** 
        \ingroup stencil_operations
    */
    inline Point ratio() const {return m_r;}
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
    
    /// Empty Query 
    /** 
        \ingroup stencil_operations
    */
    inline bool empty(){return (m_s.size() <= 0);}
    
    /// Size Query 
    /** 
        Returns the number of Stencil components in *this.
        \ingroup stencil_operations
    */
    inline int size() const {return m_s.size();}
private:
    Point m_r;                      ///< Destination refinement ratio shared by all Stencil components
    Bx m_k;                         ///< Kernel of destination shifts
    std::vector<Stencil<T>> m_s;    ///< Array of Stencils
    bool m_closed = false;          ///< Closed boolean. True if shifts / refratios of Stencils have been set (e.g. by calling close())
};

#include "StencilImplem.H"
#include "StencilDefs.H"
} //end Proto namespace
#endif
