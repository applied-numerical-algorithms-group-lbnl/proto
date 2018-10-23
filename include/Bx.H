#ifndef _PROTO_BOX_H_
#define _PROTO_BOX_H_

#include "Point.H"
#include <cstdlib> //for size_t
#include <iostream>

namespace Proto {

// forward declarations
class BxIterator;

/// An interval in DIM dimensional space.
/**
    A Bx is a region in \f$ \mathbb{Z}^{DIM} \f$ specified by two corner Point objects, <code> high </code> and <code> low </code> INCLUSIVELY.\n
    Equivalently, a Bx is a collection of \f$DIM\f$ linear intervals in \f$\mathbb{Z}^{DIM}\f$ 
    \f$[high[i] ,low[i]] : i \in [0,DIM)\f$. \n
    Bxes are written using the notation <code> [low, high] </code>.
*/
class Bx
{
public:
    typedef BxIterator iterator; //Bx::iterator aliases to BxIterator

    ///////////////////////////////////////////////////////////////////////////////////////////////  
    /** @name Constructors */
    ///@{

    /// Default Constructor
    /**
        Builds the empty Bx <code>[(-1,-1,...,-1), (0,0,...,0)]</code>
    */
    inline Bx();

    ///
    /**
        get the surrounding nodes of the box in a particular direction
    */
    inline Bx getFaceBox(int idir) const;

    ///for defining stuff with just low and high ints
    Bx (const int* a_lo, const int* a_hi)
    {
        define(a_lo, a_hi);
    }

    ///
    void define(const int* a_lo, const int* a_hi)
    {
        Point lo, hi;
        for(int idir = 0; idir < DIM; idir++)
        {
            lo[idir] = a_lo[idir];
            hi[idir] = a_hi[idir];
        }
        *this = Bx(lo, hi);
    }
    /// 2 Point Constructor
    /**
        Builds the non-trivial Bx <code>[low, high]</code>

        \param a_low    low corner
        \param a_high   high corner
    */ 
    inline Bx(const Point& a_low, const Point& a_high);
    
    /// 1 Point (Origin) Constructor
    /**
        Builds the non-trivial Bx <code>[(0,...,0), a_pt]</code>

        \param a_pt   high corner
    */ 
    inline Bx(const Point& a_pt);

    /// Copy Constructor
    inline Bx(const Bx& a_box); 
    ///@}

    ///////////////////////////////////////////////////////////////////////////////////////////////  
    /** @name Static Functions */
    ///@{
    ///  Cube Function
    /**
        Creates the Bx <code>[Point::Zeros(), Point::Ones()*(a_size-1)]</code>, 
        a cube in \f$\mathbb{Z}^{DIM}\f$ of side length a_size

        \param a_size   side length
    */
    inline static Bx Cube(int a_size);
    ///@}
    
    ///////////////////////////////////////////////////////////////////////////////////////////////  
    /** @name Accessors And Queries */
    ///@{

    /// Access Low Corner
    /**
        This function returns by value; a Bx cannot be altered through the output of <code>low()</code>.
    */
    inline Point low() const {return m_low;};

    /// Access High Corner
    /**
        This function returns by value; a Bx cannot be altered through the output of <code>high()</code>.
    */
    inline Point high() const {return m_high;}; 

    /// Edge Size
    /**
        Returns the "edge length" of this Bx on a given axis
        
        /param a_dim    direction axis in <code>[0,DIM)</code>

        Example:
        @code
        //DIM=2
        using namespace Proto;
        Bx B(Point::Zeros(), Point({1,3}));
        B.print(); //prints [(0,0),(1,3)]
        std::cout << B.size(0) << std::endl; //prints 2
        std::cout << B.size(1) << std::endl; //prints 4
        @endcode
    */
    inline std::size_t size(unsigned char a_dim) const;
    
    /// Volumetric Size
    /**
        Returns the "volume" (in number of points) of this box.
        
        Example:
        @code
        //DIM=2
        using namespace Proto;
        Bx B(Point::Zeros(), Point({1,3}));
        B.print(); //prints [(0,0),(1,3)]
        std::cout << B.size() << std::endl; //prints 8 (=2*4)
        @endcode
    */
    inline std::size_t size() const 
    {
        if (m_size <= 0){return 0;}
        else {
            return m_size;
        }
    };

    /// Contains Point Query
    /**
        Checks if a_pt is inside of *this.
        Note that <code>this->low() </code> and <code> this->high() </code> are both considered inside *this
        
        Example:
        @code
        //DIM=2
        using namespace Proto;
        Bx B(Point::Zeros(), Point({1,3}));
        std::cout << B.contains(Point({0,0})) << std::endl; //prints true
        std::cout << B.contains(Point({1,1})) << std::endl; //prints true
        std::cout << B.contains(Point({1,3})) << std::endl; //prints true
        std::cout << B.contains(Point({1,4})) << std::endl; //prints false
        @endcode
    */
    inline bool contains(const Point& a_pt) const;

    /// Contains Bx Query
    /**
        Check if Bx a_rhs is a subset of *this. Returns true even if a_rhs is not a proper subset.

        Example:
        @code
        //DIM=2
        using namespace Proto;
        Bx B0(Point({0,0}),Point({3,3}));
        Bx B1(Point({1,1}),Point({2,2}));
        Bx B2(Point({1,1}),Point({4,4}));
        Bx B2(Point({-1,-1}),Point({2,2}));
        std::cout << B0.contains(B0) << std::endl; //prints true
        std::cout << B0.contains(B1) << std::endl; //prints true
        std::cout << B0.contains(B2) << std::endl; //prints false
        std::cout << B0.contains(B3) << std::endl; //prints false
        @endcode
    */
    inline bool contains(const Bx& a_rhs) const;
    
    /// Point on Boundry Query
    /**
        Check of a_p is part of the boundry of this Bx.
    */
    inline bool onBoundary(const Point& a_p) const;
    
    /// Empty Query
    /**
        Check if *this contains no points
    */
    inline bool empty() const {return (m_size < 1);}

    /// Point to Linear Index
    /**
        Returns a linear index in <code>[0,this->size())</code> associated with <code>a_pt</code>.
        Fails by assertion if <code>a_pt</code> is not inside <code>*this</code>. 
        Indices start at <code>this->low()</code> which has an index of 0. Indices increase
        fastest in dimension 0 and slowest in dimension <code>DIM-1</code>. <code>this->high()</code>
        has an index of <code>this->size()-1</code>
        
        \param a_pt     a Point inside *this
    */
    inline unsigned int index(const Point& a_pt) const;
    
    /// Access Point by Index
    /**
        Return a Point associated with <code> a_index </code>.
        This is the inverse of the <code>index(Point)</code> function. 
        
        \param a_index  an index in <code>[0,this->size())</code>
    */
    inline Point operator[](unsigned int a_index) const;
    ///@}

    ///////////////////////////////////////////////////////////////////////////////////////////////  
    /** @name Operators */
    ///@{

    ///  Intersection Operator 
    /**
        Returns a new Bx which is the intersection of <code>*this</code> and <code>*a_rightBx</code>.

        Example:
        @code
        //DIM=2
        using namespace Proto;
        Bx B0 = Bx::Cube(3);
        Bx B1 = Bx::Cube(3).shift(Point::Ones());
        Bx B2 = B0 & B1;
        B0.print(); //prints [(0,0),(2,2)]
        B1.print(); //prints [(1,1),(3,3)]
        B2.print(); //prints [(1,1),(2,2)]
        @endcode
    */
    inline Bx operator&(const Bx& a_rightBx) const;

    inline Bx operator&(Point a_pt) const;
    
    //inline void operator&=(const Bx& a_rightBx);
  
    /// Equality Operator
    /**
        Two Bxes are considered equal if they have identical (==) low() and high()
        
        \param a_rhsBx     A Bx to test equality with
    */
    inline bool operator==(const Bx& a_rhsBx) const;
    
    /// Modulus Operator
    /**
        Convenience operator for Bx::mod(...).    
    */
    inline Point operator%(const Point& a_pt) const;
    
    ///@}

    ///////////////////////////////////////////////////////////////////////////////////////////////  
    /** @name Transformations */
    ///@{

    /// Shift Transformation
    /**
        Creates a new Bx shifted in <code>a_direction</code> by <code>a_offset</code> 
        units relative to <code>*this</code>.     
        
        \param a_direction  an int in /f$[0,DIM)/f$ specifying a direction
        \param a_offset     number of points to offset. Use a negative value for a negative offset.
    
        Example:
        @code
        //DIM=2
        using namespace Proto;
        Bx B1 = Bx::Cube(2);
        Bx B2 = B1.shift(0,1);
        Bx B3 = B1.shift(1,-1);
        B1.print(); //prints [(0,0), (1,1)]         
        B2.print(); //prints [(1,0), (2,1)]         
        B3.print(); //prints [(0,-1), (1,0)]
        @endcode         
    */
    inline Bx shift(int a_direction, int a_offset) const;

    /// Point Shift Transformation
    /**
        Creates a new Bx shifted by <code>a_pt</code> relative to <code>*this</code>.
        New Bx is: <code>[low() + a_pt, high() + a_pt]</code>
        \param a_pt     shift offset
   
        Example:
        @code
        //DIM=2
        using namespace Proto;
        Bx B1(Point::Zeros(),Point::Ones());
        Bx B2 = B1.shift(Point({2,-3}));
        B1.print(); //prints [(0,0), (1,1)]
        B2.print(); //prints [(2,-3), (3,-2)]
        @endcode
    */
    inline Bx shift(const Point& a_pt) const;

    /// Anisotropic Grow Operation
    /**
        Returns the reshaped Bx: [low() - a_pt, high() + a_pt]
        \param a_pt growth offset

        Example:
        @code
        //DIM=2
        using namespace Proto;
        Bx B0(Point::Zeros(),Point::Ones()*3);
        Bx B1 = B0.grow(Point({-1,1})); //shrink in dimension 0, grow in dimension 1
        B0.print(); //prints [(0, 0), (3, 3)]
        B1.print(); //prints [(1, -1), (2, 4)]
        @endcode
    */
    inline Bx grow(const Point& a_pt) const;

    /// Isotropic Grow Operation
    /**
        Returns a new Bx which is larger in every direction by a_numpoints units.
        If a_numpoints is negative, the box will shrink. 
        \param a_numpoints  number of points to grow by in each direction

        Example:
        @code
        //DIM=2
        using namespace Proto;
        Bx B0(Point::Zeros(),Point::Ones());
        Bx B1 = B0.grow(3);
        Bx B2 = B1.grow(-2);
        B0.print(); //prints [(0, 0), (1, 1)]
        B1.print(); //prints [(-3,-3), (4, 4)]
        B2.print(); //prints [(-1,-1), (2, 2)]
        @endcode
    */
    inline Bx grow(int a_numpoints) const;

    /// Anisotropic Coarsen Operation
    /**
        Returns a new Bx coarsened in each direction according to a_pt.
        Fails if user tries to coarsen using a non-positive ratio.

        \param a_pt    Coarsening ratios

        Example:
        @code
        //DIM=2
        using namespace Proto;
        Bx B0 = Bx::Cube(4);
        //coarsen in only the 1 direction
        Bx B1 = B0.coarsen(Point({1,2}));
        B1.print(); //prints [(0,0),(3,1)]
        @endcode
    */
    inline Bx coarsen(const Point& a_pt) const;

    /// Isotropic Coarsen Operation
    /**
        Returns a new Bx coarsened by a factor of a_numpoints. Fails if user tries to 
        coarsen using a non-positive ratio.

        \param a_ratio     Coarsening ratio

        Example:
        @code
        //DIM=2
        using namespace Proto;
        Bx B0 = Bx::Cube(4);
        Bx B1 = B0.coarsen(2);
        B0.print(); //prints [(0, 0), (3, 3)]
        B1.print(); //prints [(0, 0), (1, 1)]

        Bx B2 = Bx::Cube(3).shift(Point::Ones()*2);
        Bx B3 = B2.coarsen(2);
        B2.print(); //prints [(2, 2), (4, 4)]
        B3.print(); //prints [(1, 1), (2, 2)]
        @endcode
    */
    inline Bx coarsen(unsigned int a_ratio) const;

    /// return true if the box is coarsenable by the input integer.
    inline bool coarsenable(const int& a_ratio) const
    {
        Bx coarBox  =  this->coarsen(Point::Ones(a_ratio));
        Bx refCoBox = coarBox.refine(Point::Ones(a_ratio));
        return(refCoBox == (*this));
    }

    // Same as coarsen, but the lower bound is rounded UP. Used in Stencil
    inline Bx taperCoarsen(const Point& a_ref)
    {
        Point H = m_high/a_ref;
        Point L;
        for (int ii = 0; ii < DIM; ii++)
        {
            L[ii] = ceil(m_low[ii]/((double)a_ref[ii]));
        }
        return Bx(L,H);
    }

    /// Anisotropic Refine Operation
    /**
        Returns a new Bx refined in each direction according to a_pt.
        \param a_pt    Refinement ratios

        Example:
        @code
        //DIM=2
        using namespace Proto;
        Bx B0 = Bx::Cube(2);
        Bx B1 = B0.refine(Point({1,2}));
        B0.print(); //prints [(0, 0), (1, 1)]
        B1.print(); //prints [(0, 0), (1, 3)]
        @endcode
        
    */
    inline Bx refine(const Point& a_pt) const;
    
    /// Isotropic Refine Operation
    /**
        Returns a new Bx refined in all direcitons by a_numpoints.
        \param a_ratio     Refinement ratio

        Example:
        @code
        //DIM=2
        using namespace Proto;
        Bx B0 = Bx::Cube(2);
        Bx B1 = B0.refine(2);
        B0.print(); //prints [(0, 0), (1, 1)]
        B1.print(); //prints [(0, 0), (3, 3)]

        Bx B2 = Bx::Cube(2).shift(Point::Ones());
        Bx B3 = B2.coarsen(2);
        B2.print(); //prints [(1, 1), (2, 2)]
        B3.print(); //prints [(2, 2), (5, 5)]
        @endcode
    */
    inline Bx refine(unsigned int a_ratio) const;

    /// Flatten
    /**
        Returns a copy of *this with dimension a_dir flattened to a thickness of 1.
        Useful for creating Bxes of dimensionality less than DIM (e.g. a plane in 3D)

        \param a_dir    Direction to flatten
        \param a_upper  Flatten upwards?
        Example:
        @code
        //DIM=3
        using namespace Proto;
        Bx B0 = Bx::Cube(4).shift(Point::Ones());
        Bx B1 = B0.flatten(2);
        Bx B2 = B0.flatten(2,true);
        B0.print(); //prints [(1, 1, 1), (4, 4, 4)]
        B1.print(); //prints [(1, 1, 1), (4, 4, 1)]
        B2.print(); //prints [(1, 1, 4), (4, 4, 4)]
        @endcode
    */
    inline Bx flatten(const int a_dir, bool a_upper = false) const;

    
    /// Edge
    /**
        Returns the subset on the boundary of *this in a given direction.
        Very handy for adjacency and boundary based computations. 

        this->edge(Point::Basis(0,+1)) returns the upper face along the 0 axis.
        this->edge(Point::Basis(0,-1)) returns the lower face along the 0 axis.
        this->edge(Point::Ones()) returns m_high
        this->edge(Point::Ones(-1) returns m_low
        this->edge(Point(1,-1,0,...,0) returns the 1D line of points which forms the intersection
          of the upper face along the 0 axis and the lower face along the 1 axis. 
    */
    inline Bx edge(const Point& a_dir) const;
    
    /// Extrude
    /**
        Returns a copy of *this with dimension a_dir extruded a distance a_dist.
        By default, the extrusion is downwards. If a_up is true, extrudes upwards.  
        
        \param a_dir    Direction to extrude
        \param a_dist   Distance to extrude
        \param a_upper  Extrude upwards?
        Example:
        @code
        //DIM=3
        using namespace Proto;
        Bx B0 = Bx::Cube(4).shift(Point::Ones());
        Bx B1 = B0.extrude(2,3);
        Bx B2 = B0.extrude(2,3,true);
        B0.print(); //prints [(1, 1, 1), (4, 4, 4)]
        B1.print(); //prints [(1, 1, -2), (4, 4, 4)]
        B2.print(); //prints [(1, 1, 1), (4, 4, 7)]
        @endcode

    */
    inline Bx extrude(int a_dir, int a_dist, bool a_upper = false) const;
    ///@}
  
    
    ///////////////////////////////////////////////////////////////////////////////////////////////  
    /** @name Utility */
    ///@{
    
    /// Iterator Begin
    /**
        See documentation for Proto::BxIterator for a basic usage example.
    */
    inline BxIterator begin() const; 
    
    /// Iterator End
    /**
        See documentation for Proto::BxIterator for a basic usage example.
    */
    inline BxIterator end() const; 
    
    /// Iterator Reverse Begin
    /**
        See documentation for Proto::BxIterator for a basic usage example.
    */
    inline BxIterator rbegin() const; 
    
    /// Iterator Reverse End
    /**
        See documentation for Proto::BxIterator for a basic usage example.
    */
    inline BxIterator rend() const; 
    
    /// Point Modulus Operation
    /**
        Returns the periodic image of a_pt that is inside *this.
        \param a_pt   Point divisor

        Example:
        @code
        //DIM = 2;
        using namespace Proto;
        Bx B = Bx::Cube(3);
        std::cout << B.mod(Point({3,3})) << std::endl; //prints (1,1)
        std::cout << B.mod(Point({-1,-1})) << std::endl; //prints (2,2)
        @endcode
    */
    inline Point mod(const Point& a_pt) const;

    /// Print
    /**
        Prints *this using the format <code>[low, high]
    */ 
    inline void print() const; 
    ///@}

private:  
    inline void recomputeSize(); ///< Used to reevaluate the size of the box when it is changed.
    Point m_low; ///< Point object containing the lower bounds of the Bx.
    Point m_high;  ///< Point object containing the upper bounds of the Bx.
    int m_size; ///< "Volume" of the box.
}; //end class Bx

/// OStream Operator 
inline std::ostream& operator<<(std::ostream& a_os, const Bx& a_box)
{
    a_os << "[" << a_box.low() << ", " << a_box.high() << "]";
    return a_os;
}

/// Iterator for Bxes
/**
    Iteration class which conforms to most of the syntax of std Iterator implementations.

    Example:
    @code
    using namespace Proto;
    int domainSize = 16;
    Bx B = Bx::Cube(domainSize);
    for (auto iter = B.begin(); iter != B.end(); ++iter)
    {
        std::cout << *iter << std::endl; //prints current Point
        std::cout << iter() << std::endl; //prints current index in [0,B.size())
        
        //... do stuff with *iter...
    }
    @endcode
*/
class BxIterator
{
public:
    /// Constructor
    
    BxIterator(const Bx& a_box, int a_pos = 0)
    : m_box(a_box), m_pos(a_pos) {}
    
    BxIterator(const Bx& a_box, const Point& a_pos)
    : m_box(a_box)
    {
        m_pos = m_box.index(a_pos);
    }

    /// Equality Operator
    /**
        Determines equality between two BxIterators. Note that this operation CAN return true
        for BxIterators constructed from different Bx instances so long as said Bxes have the
        same bounds (e.g. they must be ==)
    */
    inline bool operator==(const BxIterator& a_iter) const
    {
        return ((m_pos == a_iter.m_pos) && (m_box == a_iter.m_box));
    };

    /// Inequality Operator
    inline bool operator!=(const BxIterator& a_iter) const
    {
        return !((*this) == a_iter);
    };
    
    /// Dereference Iterator
    /**
        Returns a Point
    */
    inline Point operator*() const
    {
        return m_box[m_pos];
    };

    inline int operator()() const
    {
        return m_pos;
    };

    /// Increment Iterator
    inline BxIterator& operator++()
    {
        m_pos++;
        return *this;
    };
    
    /// Decrement Iterator
    inline BxIterator& operator--()
    {
        m_pos--;
        return *this;
    };

    /// Check if Done
    inline bool done()
    {
        return ((m_pos < 0) && (m_pos >= m_box.size()));
    }

private:
    const Bx m_box;
    int m_pos;
}; //end class BxIterator

#include "BxImplem.H"
} //end namespace Proto
#endif
