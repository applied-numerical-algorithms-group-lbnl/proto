#ifndef _PROTO_BOX_H___
#define _PROTO_BOX_H___

#include "Proto_Point.H"
#include <cstdlib> //for size_t
#include <iostream>
#include <cmath>
using Proto::Point;
namespace Proto 
{

// forward declarations
  class BoxIterator;

/// An interval in DIM dimensional space.
/**
   A Box is a region in \f$ \bf{Z}^{DIM} \f$ specified by two corner Point objects, <code> high </code> and <code> low </code> INCLUSIVELY.\n
   Equivalently, a Box is a collection of \f$DIM\f$ linear intervals in \f$\bf{Z}^{DIM}\f$ 
   \f$[high[i] ,low[i]] : i \in [0,DIM)\f$. \n
   Boxes are written using the notation <code> [low, high] </code>.
*/
  class Box
  {
  public:
    typedef BoxIterator iterator; //Box::iterator aliases to BoxIterator

    ///////////////////////////////////////////////////////////////////////////////////////////////  
    /** @name Constructors */
    ///@{

    /// Default Constructor
    /**
       Builds the empty Box <code>[(-1,-1,...,-1), (0,0,...,0)]</code>
    */
    inline Box();
    
    /// 2 Point Constructor
    /**
       Builds the non-trivial Bx <code>[low, high]</code>

       \param a_low    low corner
       \param a_high   high corner
    */ 
    inline Box(const Point& a_low, const Point& a_high);

    /// Primitive Constructor
    /**
       Used for building a Bx on platforms where plain-old-data is more convenient
       \param a_lo   C-Array representing this.low()
       \param a_hi   C-Array representing this.high()
    */
    inline Box (const int* a_lo, const int* a_hi)
    {
      define(a_lo, a_hi);
    }

    /// Define
    /**
       Weak construct this Bx if necessary
      
       \param a_lo   C-Array representing this.low()
       \param a_hi   C-Array representing this.high()
    */
    void define(const int* a_lo, const int* a_hi)
    {
      Point lo, hi;
      for(int idir = 0; idir < DIM; idir++)
      {
        lo[idir] = a_lo[idir];
        hi[idir] = a_hi[idir];
      }
      *this = Box(lo, hi);
    }
    
    /// 1 Point (Origin) Constructor
    /**
       Builds the non-trivial Box <code>[(0,...,0), a_pt]</code>

       \param a_pt   high corner
    */ 
    inline Box(const Point& a_pt);

    /// Copy Constructor
    inline Box(const Box& a_box); 
    ///@}

    ///////////////////////////////////////////////////////////////////////////////////////////////  
    /** @name Static Functions */
    ///@{
    ///  Cube Function
    /**
       Creates the Box <code>[Point::Zeros(), Point::Ones()*(a_size-1)]</code>, 
       a cube in \f$\bf{Z}^{DIM}\f$ of side length a_size

       \param a_size   side length
    */
    inline static Box Cube(int a_size);
    
    ///  Kernel Function
    /**
        Creates a Box of size <code> 2*a_radius + 1</code> on all sides with (0,...,0) at the center.
        Useful for iterating through boundary cells.

        \param a_radius the number of "rings" of cells around the center cell. 
    */
    inline static Box Kernel(int a_radius);
    ///@}
    
    ///////////////////////////////////////////////////////////////////////////////////////////////  
    /** @name Accessors And Queries */
    ///@{

    /// Access Low Corner
    /**
       This function returns by value; a Box cannot be altered through the output of <code>low()</code>.
    */
    inline Point low() const {return m_low;};

    /// Access High Corner
    /**
       This function returns by value; a Box cannot be altered through the output of <code>high()</code>.
    */
    inline Point high() const {return m_high;}; 

    /// Edge Size
    /**
       Returns the "edge length" of this Box on a given axis
        
       /param a_dim    direction axis in <code>[0,DIM)</code>

       Example:
       @code
       //DIM=2
       using namespace Proto;
       Box B(Point::Zeros(), Point({1,3}));
       B.print(); //prints [(0,0),(1,3)]
       std::cout << B.size(0) << std::endl; //prints 2
       std::cout << B.size(1) << std::endl; //prints 4
       @endcode
    */
    CUDA_DECORATION
    inline std::size_t size(unsigned char a_dim) const;
    
    /// Volumetric Size
    /**
       Returns the "volume" (in number of points) of this box.
        
       Example:
       @code
       //DIM=2
       using namespace Proto;
       Box B(Point::Zeros(), Point({1,3}));
       B.print(); //prints [(0,0),(1,3)]
       std::cout << B.size() << std::endl; //prints 8 (=2*4)
       @endcode
    */
    CUDA_DECORATION
    inline std::size_t size() const 
    {
      if (m_size <= 0)
      {
        return 0;
      }
      else 
      {
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
       Box B(Point::Zeros(), Point({1,3}));
       std::cout << B.contains(Point({0,0})) << std::endl; //prints true
       std::cout << B.contains(Point({1,1})) << std::endl; //prints true
       std::cout << B.contains(Point({1,3})) << std::endl; //prints true
       std::cout << B.contains(Point({1,4})) << std::endl; //prints false
       @endcode
    */
    inline bool contains(const Point& a_pt) const;

    /// Contains Box Query
    /**
       Check if Bx a_rhs is a subset of *this. Returns true even if a_rhs is not a proper subset.
       Equivalent to <code> this->contains(a_rhs.low()) && this->contains(a_rhs.high() </code>
        
       Example:
       @code
       //DIM=2
       using namespace Proto;
       Box B0(Point({0,0}),Point({3,3}));
       Box B1(Point({1,1}),Point({2,2}));
       Box B2(Point({1,1}),Point({4,4}));
       Box B2(Point({-1,-1}),Point({2,2}));
       std::cout << B0.contains(B0) << std::endl; //prints true
       std::cout << B0.contains(B1) << std::endl; //prints true
       std::cout << B0.contains(B2) << std::endl; //prints false
       std::cout << B0.contains(B3) << std::endl; //prints false
       @endcode
    */
    inline bool contains(const Box& a_rhs) const;
    
    /// Point on Boundry Query
    /**
       Check of a_p is part of the boundry of this Bx.

       \param a_p  A Point
    */
    inline bool onBoundary(const Point& a_p) const;
    
    /// Empty Query
    /**
       Check if *this contains no Points
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
    CUDA_DECORATION
    inline unsigned int index(const Point& a_pt) const;
    
    /// Access Point by Index
    /**
       Return a Point associated with <code> a_index </code>.
       This is the inverse of the <code>index(Point)</code> function. 
        
       \param a_index  an index in <code>[0,this->size())</code>
    */
    CUDA_DECORATION
    inline Point operator[](unsigned int a_index) const;

    ///@}


    ///for stencil internals
    CUDA_DECORATION
    static unsigned int sindex(Point a_pt, Point a_low, Point a_high) 
    {
      unsigned int factor = 1;
      unsigned int linIndex = a_pt[0] - a_low[0];
      for (unsigned char i = 1;i < DIM;i++)
      { 
        factor = factor*(a_high[i-1] - a_low[i-1] + 1);
        linIndex = linIndex + (a_pt[i] - a_low[i])*factor;
      }
      return linIndex;
    }

    CUDA_DECORATION
    static int soffset(const Point& a_pt, const Point& a_low, const Point& a_high)
    {
      int factor = 1;
      int linIndex = a_pt[0];
      for(unsigned char i=1; i<DIM; i++)
        {
          factor = factor*(a_high[i-1] - a_low[i-1] + 1);
          linIndex = linIndex + (a_pt[i])*factor;
        }
      return linIndex;
    }     
    ///////////////////////////////////////////////////////////////////////////////////////////////  
    /** @name Operators */
    ///@{

    ///  Intersection Operator 
    /**
       Returns a new Box which is the intersection of <code>*this</code> and <code>*a_rightBox</code>.

       \param a_rightBx  Another Bx

       Example:
       @code
       //DIM=2
       using namespace Proto;
       Box B0 = Box::Cube(3);
       Box B1 = Box::Cube(3).shift(Point::Ones());
       Box B2 = B0 & B1;
       B0.print(); //prints [(0,0),(2,2)]
       B1.print(); //prints [(1,1),(3,3)]
       B2.print(); //prints [(1,1),(2,2)]
       @endcode
    */
    inline Box operator&(const Box& a_rightBox) const;

    inline void operator&=(const Box& a_rightBox) 
    {
      Box interbox = a_rightBox & (*this);
      *this = interbox;
    }

    // TODO: Probably want to change this syntax to an actual function; it's potentially confusing with intersection
    /// Add Point Operator
    /**
       Returns a new Bx which has been grown such that it contains a_pt. 

       \param a_pt   Point to be added
    */
    inline Box operator&(Point a_pt) const;
    
    //inline void operator&=(const Box& a_rightBox);
  
    /// Equality Operator
    /**
       Two Boxes are considered equal if they have identical (==) low() and high()
        
       \param a_rhsBox     A Box to test equality with
    */
    inline bool operator==(const Box& a_rhsBox) const;

    /// Inequality Operator
    inline bool operator!=(const Box& a_rhsBox) const
    {
      bool equals = (*this==a_rhsBox);
      return !equals;
    }
    
    /// Modulus Operator
    /**
       Convenience operator for Bx::mod(...).

       \param a_pt   A Point to mod by    
    */
    inline Point operator%(const Point& a_pt) const;
    
    ///@}

    ///////////////////////////////////////////////////////////////////////////////////////////////  
    /** @name Transformations */
    ///@{

    /// Shift Transformation
    /**
       Creates a new Box shifted in <code>a_direction</code> by <code>a_offset</code> 
       units relative to <code>*this</code>.     
        
       \param a_direction  an int in /f$[0,DIM)/f$ specifying a direction
       \param a_offset     number of points to offset. Use a negative value for a negative offset.
    
       Example:
       @code
       //DIM=2
       using namespace Proto;
       Box B1 = Box::Cube(2);
       Box B2 = B1.shift(0,1);
       Box B3 = B1.shift(1,-1);
       B1.print(); //prints [(0,0), (1,1)]         
       B2.print(); //prints [(1,0), (2,1)]         
       B3.print(); //prints [(0,-1), (1,0)]
       @endcode         
    */
    inline Box shift(int a_direction, int a_offset) const;

    /// Point Shift Transformation
    /**
       Creates a new Box shifted by <code>a_pt</code> relative to <code>*this</code>.
       New Box is: <code>[low() + a_pt, high() + a_pt]</code>
       \param a_pt     shift offset
   
       Example:
       @code
       //DIM=2
       using namespace Proto;
       Box B1(Point::Zeros(),Point::Ones());
       Box B2 = B1.shift(Point({2,-3}));
       B1.print(); //prints [(0,0), (1,1)]
       B2.print(); //prints [(2,-3), (3,-2)]
       @endcode
    */
    inline Box shift(const Point& a_pt) const;

    /// Anisotropic Grow Operation
    /**
       Returns the reshaped Box: [low() - a_pt, high() + a_pt]
       \param a_pt growth offset

       Example:
       @code
       //DIM=2
       using namespace Proto;
       Box B0(Point::Zeros(),Point::Ones()*3);
       Box B1 = B0.grow(Point({-1,1})); //shrink in dimension 0, grow in dimension 1
       B0.print(); //prints [(0, 0), (3, 3)]
       B1.print(); //prints [(1, -1), (2, 4)]
       @endcode
    */
    inline Box grow(const Point& a_pt) const;

    /// Anisotropic Grow (Int Version)
    /**
        Returns grow(Point::Basis(dir, offset)
    */
    inline Box grow(int a_dir, int a_offset) const
    {
        return grow(Point::Basis(a_dir, a_offset));
        /*
      if(igrow > 0)
      {
        return grow(Point::Basis(idir, igrow));
      }
      else
      {
        return grow(-Point::Basis(idir, igrow));
      }
      */
    }

    /// grows only high side 
    inline Box growHi(int idir, int igrow) const
    {
      Box retval(m_low, m_high + Point::Basis(idir, igrow));
      return retval;
    }
  

    ///grows only low side
    inline Box growLo(int idir, int igrow) const
    {
      Box retval(m_low - Point::Basis(idir, igrow), m_high);
      return retval;
    }
  

    /// Isotropic Grow Operation
    /**
       Returns a new Box which is larger in every direction by a_numpoints units.
       If a_numpoints is negative, the box will shrink. 
       \param a_numpoints  number of points to grow by in each direction

       Example:
       @code
       //DIM=2
       using namespace Proto;
       Box B0(Point::Zeros(),Point::Ones());
       Box B1 = B0.grow(3);
       Box B2 = B1.grow(-2);
       B0.print(); //prints [(0, 0), (1, 1)]
       B1.print(); //prints [(-3,-3), (4, 4)]
       B2.print(); //prints [(-1,-1), (2, 2)]
       @endcode
    */
    inline Box grow(int a_numpoints) const;

    /// Anisotropic Coarsen Operation
    /**
       Returns a new Box coarsened in each direction according to a_pt.
       Fails if user tries to coarsen using a non-positive ratio.

       \param a_pt    Coarsening ratios

       Example:
       @code
       //DIM=2
       using namespace Proto;
       Box B0 = Box::Cube(4);
       //coarsen in only the 1 direction
       Box B1 = B0.coarsen(Point({1,2}));
       B1.print(); //prints [(0,0),(3,1)]
       @endcode
    */
    inline Box coarsen(const Point& a_pt) const;

    /// Isotropic Coarsen Operation
    /**
       Returns a new Box coarsened by a factor of a_numpoints. Fails if user tries to 
       coarsen using a non-positive ratio.
       If the limits of *this are not multiples of the coarsening ratio (e.g. if this-> coarsenable() != true)
       the resulting Bx is not guaranteed to be a subset of *this. In this situation, use Bx::taperCoarsen

       \param a_ratio     Coarsening ratio

       Example:
       @code
       //DIM=2
       using namespace Proto;
       Box B0 = Box::Cube(4);
       Box B1 = B0.coarsen(2);
       B0.print(); //prints [(0, 0), (3, 3)]
       B1.print(); //prints [(0, 0), (1, 1)]

       Box B2 = Box::Cube(3).shift(Point::Ones()*2);
       Box B3 = B2.coarsen(2);
       B2.print(); //prints [(2, 2), (4, 4)]
       B3.print(); //prints [(1, 1), (2, 2)]
       @endcode
    */
    inline Box coarsen(unsigned int a_ratio) const;

    /// Coarsenable Query
    /**
       Returns true if *this will coarsen normally by a given ratio. 
    
       \param a_ratio    Coarsening ratio
    */
    inline bool coarsenable(const int& a_ratio) const
    {
      Box coarBox  =  this->coarsen(Point::Ones(a_ratio));
      Box refCoBox = coarBox.refine(Point::Ones(a_ratio));
      return(refCoBox == (*this));
    }

    // Tapered Coarsen
    /**
       This function is identical to Bx::coarsen when Bx::coarsenable() is true.
       For non-coarsenable Boxes, the new Box's limits are rounded such that the result
       Is always a subset of *this. Specifically, the lower limits are always rounded UP.

       \param a_ref  Coarsening ratio
    */
    inline Box taperCoarsen(const Point& a_ref)
    {
      Point H = m_high/a_ref;
      Point L;
      for (int ii = 0; ii < DIM; ii++)
      {
        L[ii] = ceil(m_low[ii]/((double)a_ref[ii]));
      }
      return Box(L,H);
    }

    /// Anisotropic Refine Operation
    /**
       Returns a new Box refined in each direction according to a_pt.
       \param a_pt    Refinement ratios

       Example:
       @code
       //DIM=2
       using namespace Proto;
       Box B0 = Box::Cube(2);
       Box B1 = B0.refine(Point({1,2}));
       B0.print(); //prints [(0, 0), (1, 1)]
       B1.print(); //prints [(0, 0), (1, 3)]
       @endcode
        
    */
    inline Box refine(const Point& a_pt) const;
    
    /// Isotropic Refine Operation
    /**
       Returns a new Box refined in all direcitons by a_numpoints.
       \param a_ratio     Refinement ratio

       Example:
       @code
       //DIM=2
       using namespace Proto;
       Box B0 = Box::Cube(2);
       Box B1 = B0.refine(2);
       B0.print(); //prints [(0, 0), (1, 1)]
       B1.print(); //prints [(0, 0), (3, 3)]

       Box B2 = Box::Cube(2).shift(Point::Ones());
       Box B3 = B2.coarsen(2);
       B2.print(); //prints [(1, 1), (2, 2)]
       B3.print(); //prints [(2, 2), (5, 5)]
       @endcode
    */
    inline Box refine(unsigned int a_ratio) const;

    /// Flatten
    /**
       Returns a copy of *this with dimension a_dir flattened to a thickness of 1.
       Useful for creating Boxes of dimensionality less than DIM (e.g. a plane in 3D)

       \param a_dir    Direction to flatten
       \param a_upper  Flatten upwards?
       Example:
       @code
       //DIM=3
       using namespace Proto;
       Box B0 = Box::Cube(4).shift(Point::Ones());
       Box B1 = B0.flatten(2);
       Box B2 = B0.flatten(2,true);
       B0.print(); //prints [(1, 1, 1), (4, 4, 4)]
       B1.print(); //prints [(1, 1, 1), (4, 4, 1)]
       B2.print(); //prints [(1, 1, 4), (4, 4, 4)]
       @endcode
    */
    inline Box flatten(const int a_dir, bool a_upper = false) const;


    ///
    inline Box adjCellLo(int a_idir, int a_length) const
    {
      Point lo = low();
      Point hi = high();
      int sm = lo[a_idir];

      lo[a_idir] = sm - a_length;
      hi[a_idir] = sm - 1;
      return Box(lo, hi);
    }

    ///
    inline Box adjCellHi(int a_idir, int a_length) const
    {
      Point lo = low();
      Point hi = high();
      int bg = hi[a_idir]+1;

      lo[a_idir] = bg;
      hi[a_idir] = bg + a_length - 1;
      return Box(lo, hi);
    }

    /// Edge
    /**
       Returns the subset on the boundary of *this in a given direction.
       Very handy for adjacency and boundary based computations. 
       
       \param a_dir  Direction of desired edge with respect to the cell center of *this
 
       Examples:
       @code
       this->edge(Point::Basis(0,+1)); //returns the upper face along the 0 axis.
       this->edge(Point::Basis(0,-1)); //returns the lower face along the 0 axis.
       this->edge(Point::Ones());      //returns m_high
       this->edge(Point::Ones(-1);     //returns m_low
       this->edge(Point(1,-1,0,...,0); //returns the 1D line of points which forms the intersection
       //of the upper face along the 0 axis and the lower face along the 1 axis. 
       @endcode
    */
    inline Box edge(const Point& a_dir, int a_dist = 1) const;
   
    /// Adjacent
    /**
        Returns a box adjacent to *this in a given direction with a thickness
        a_dist in the normal direction. If there are multiple non-zero entries
        in a_dir, a_dist will be applied to the thickness in all of those directions

        \param a_dir "Normal" direction of the desired adjacent cell
        \param a_dist "thickness" of the desired adjacent cell

        Examples:
        @code
        //DIM = 2;
        Box B0 = Box::Cube(8);  // [(0,0) ,  (7,7)]
        B0.adjacent((1,0) , 2); // [(8,0) ,  (9,7)]
        B0.adjacent((0,-1), 2); // [(0,-2), (7,-1)]
        B0.adjacent((-1,1), 2);  // [(-2,8), (-1,9)]
        @endcode 
    */
    inline Box adjacent(const Point& a_dir, int a_dist = -1) const;
    
    /// Extrude
    /**
       Returns a new Bx with the <code> a_dir </code> direction extruded a distance <code> a_dist </code>.
       By default, the extrusion is upwards. If a_up is false, extrudes downwards. Useful computing e.g. "face" boxes.
        
       \param a_dir    Axis to extrude along
       \param a_dist   (Optional) Distance to extrude (default: 1)
       \param a_upper  (Optional) Extrude upwards? (default: true)

       Example:
       @code
       //DIM=3
       using namespace Proto;
       Box B0 = Box::Cube(4).shift(Point::Ones());
       Box B1 = B0.extrude(2,3);
       Box B2 = B0.extrude(2,3,true);
       B0.print(); //prints [(1, 1, 1), (4, 4, 4)]
       B1.print(); //prints [(1, 1, -2), (4, 4, 4)]
       B2.print(); //prints [(1, 1, 1), (4, 4, 7)]
       @endcode

    */
    inline Box extrude(int a_dir, int a_dist = 1, bool a_upper = true) const;
    ///@}
  
    ///////////////////////////////////////////////////////////////////////////////////////////////  
    /** @name Utility */
    ///@{
    
    /// Iterator Begin
    /**
       See documentation for Proto::BoxIterator for a basic usage example.
    */
    inline BoxIterator begin() const; 
    
    /// Iterator End
    /**
       See documentation for Proto::BoxIterator for a basic usage example.
    */
    inline BoxIterator end() const; 
    
    /// Iterator Reverse Begin
    /**
       See documentation for Proto::BoxIterator for a basic usage example.
    */
    inline BoxIterator rbegin() const; 
    
    /// Iterator Reverse End
    /**
       See documentation for Proto::BoxIterator for a basic usage example.
    */
    inline BoxIterator rend() const; 
    
    /// Point Modulus Operation
    /**
       Returns the periodic image of a_pt that is inside *this.
       \param a_pt   Point divisor

       Example:
       @code
       //DIM = 2;
       using namespace Proto;
       Box B = Box::Cube(3);
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
    Point m_low; ///< Point object containing the lower bounds of the Box.
    Point m_high;  ///< Point object containing the upper bounds of the Box.
    int m_size; ///< "Volume" of the box.
  }; //end class Box

/// OStream Operator 
  inline std::ostream& operator<<(std::ostream& a_os, const Box& a_box)
  {
    a_os << "[" << a_box.low() << ", " << a_box.high() << "]";
    return a_os;
  }

/// Iterator for Boxes
/**
   Iteration class which conforms to most of the syntax of std Iterator implementations.

   Example:
   \snippet Snippets.cpp proto_boxiter
*/
  class BoxIterator
  {
  public:
    /// Constructor
    
    BoxIterator(const Box& a_box, int a_pos = 0)
      : m_box(a_box), m_pos(a_pos) {}
    
    BoxIterator(const Box& a_box, const Point& a_pos)
      : m_box(a_box)
    {
      m_pos = m_box.index(a_pos);
    }

    /// Equality Operator
    /**
       Determines equality between two BoxIterators. Note that this operation CAN return true
       for BoxIterators constructed from different Box instances so long as said Boxes have the
       same bounds (e.g. they must be ==)
    */
    inline bool operator==(const BoxIterator& a_iter) const
    {
      return ((m_pos == a_iter.m_pos) && (m_box == a_iter.m_box));
    };

    /// Inequality Operator
    inline bool operator!=(const BoxIterator& a_iter) const
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

    /// Get Index
    /**
       Returns the current index of the point *(*this).  
    */
    inline int operator()() const
    {
      return m_pos;
    };

    /// Increment Iterator
    inline BoxIterator& operator++()
    {
      m_pos++;
      return *this;
    };
    
    /// Decrement Iterator
    inline BoxIterator& operator--()
    {
      m_pos--;
      return *this;
    };

    /// Check if Done
    inline bool done()
    {
      return ((m_pos < 0) || (m_pos >= m_box.size()));
    }


  private:
    const Box m_box;
    int m_pos;
  }; //end class BoxIterator

#include "implem/Proto_BoxImplem.H"
} //end namespace Proto
#endif
