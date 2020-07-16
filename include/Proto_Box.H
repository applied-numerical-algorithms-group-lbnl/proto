#ifndef _PROTO_BOX_H___
#define _PROTO_BOX_H___

#include "Proto_Point.H"
#include "Proto_LoHiSide.H"
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

    //useful for bc evaluation
    Point boundary(const Side::LoHiSide a_side) const
    {
      Point retval = low();
      if(a_side == Side::Hi)
      {
        retval =  high();
      }
      return retval;
    }

    //useful for bc evaluation
    Box faceBox(int a_idir, const Side::LoHiSide a_side) const
    {
      Box retval = adjCellSide(a_idir, 1, a_side);
      if(a_side == Side::Lo)
      {
        retval =  retval.shift(Point::Basis(a_idir, 1));
      }
      return retval;
    }
    
    /// 2 Point Constructor
    /**
       Builds the non-trivial Box <code>[low, high]</code>

       \param a_low    low corner
       \param a_high   high corner
    */ 
    inline Box(const Point& a_low, const Point& a_high);

    /// Primitive Constructor
    /**
       Used for building a Box on platforms where plain-old-data is more convenient
       \param a_lo   C-Array representing this.low()
       \param a_hi   C-Array representing this.high()
    */
    inline Box (const int* a_lo, const int* a_hi) { define(a_lo, a_hi); }

    /// Define
    /**
       Weak construct this Box if necessary
      
       \param a_lo   C-Array representing this.low()
       \param a_hi   C-Array representing this.high()
    */
    inline void define(const int* a_lo, const int* a_hi);
    
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
        
        Example:
        @code
        //DIM=2
        using namespace Proto;
        auto K0 = Box::Kernel(1);
        K0.print(); // [(-1,-1), (1, 1)] 
        auto K1 = Box::Kernel(2);
        K1.print(); // [(-2, -2), (2, 2)]
        @endcode
    */
    inline static Box Kernel(int a_radius);
   
    /// Stencil Index
    /**
        @private
        Returns the linear index of <code>a_pt</code> inside the <code>Box</code> defined by <code>[a_low, a_high]. Equivalent to <code>Box::index</code> but does not require a <code>Box</code> instance.
    
        \param a_pt     A Point
        \param a_low    The low Point of a Box
        \param a_high   The high Point of a Box
    */
    CUDA_DECORATION
    inline static unsigned int sindex(Point a_pt, Point a_low, Point a_high); 
   
    /// Stencil Offset
    /**
        @private
        Used internally by <code>Stencil</code>, not recommended for public use.
        
        \param a_pt     A Point
        \param a_low    The low Point of a Box
        \param a_high   The high Point of a Box
    */
    CUDA_DECORATION
    inline static unsigned int soffset(const Point& a_pt, const Point& a_low, const Point& a_high);
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
    inline std::size_t size() const; 

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
       Check if Box a_rhs is a subset of *this. Returns true even if a_rhs is not a proper subset.
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
       Check of a_p is part of the boundry of this Box.

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



    ///////////////////////////////////////////////////////////////////////////////////////////////  
    /** @name Operators */
    ///@{

    ///  Intersection Operator 
    /**
       Returns a new Box which is the intersection of <code>*this</code> and <code>*a_rightBox</code>.

       \param a_rightBox  Another Box

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

    
    /// In Place Intersection
    /**
        Performs intersection on *this with <code>a_rhs</code>
        
        \param a_rhs    Another Box 
    */
    inline void operator&=(const Box& a_rhs); 

    /// Add Point Operator
    /**
        Returns a new Box which has been grown such that it contains a_pt. 
        if <code>this->contains(a_pt></code> this function returns *this.

        \param a_pt   Point to be added
    
        Example:
        @code
        //DIM=2
        using namespace Proto;
        
        auto B0 = Box::Cube(3);
        B0.print(); //[(0,0), (2,2)]
        Point p(-1,3);
        auto B1 = B0 & p;
        B1.print(); //[(-1,0), (2,3)]
        @endcode
    */
    inline Box operator&(const Point& a_pt) const;
  
    /// In Place Add Point Operator
    /**
        Alters <code>*this</code> so that it now includes <code>a_pt</code>.
        If <code>this->contains(a_pt)</code>, this function does nothing.

        \param a_pt     A Point
    */
    inline void operator&=(const Point& a_pt);
    
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
       Convenience operator for Box::mod(...).
       Returns the periodic image of <code>a_pt</code> contained in <code>*this</code>

       \param a_pt   A Point to mod by    
    */
    inline Point operator%(const Point& a_pt) const;
    
    /// Modulus Function
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
    
    /// Anisotropic Grow Operation
    /**
       Returns the reshaped Box: [low() - a_pt, high() + a_pt]
       To grow a grow a box only upwards / downwards, see <code>Box::extrude</code>.
       
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
        Returns grow(Point::Basis(dir, offset))
        Useful in places where plain-old-data inputs are preferable.
       To grow a grow a box only upwards / downwards, see <code>Box::extrude</code>.

        \param a_dir     A direction in [0,DIM)
        \param a_offset  Distance to grow
    */
    inline Box grow(int a_dir, int a_offset) const;
    
    /// Extrude
    /**
       AKA single sided grow. Returns a new Box with the <code> a_dir </code> direction extruded a distance <code> a_dist </code).
       By default, the extrusion is upwards. If <code>a_upper</code> is false, extrudes downwards.
        
       \param a_dir    Axis to extrude along
       \param a_dist   (Optional) Distance to extrude (default: 1)

       Example:
       @code
       //DIM=3
       using namespace Proto;
       Box B0 = Box::Cube(4).shift(Point::Ones());
       Box B1 = B0.extrude(Point::Ones(), 2);
       Box B2 = B0.extrude(Point::Basis(0,-1),3);
       B0.print(); //prints [(1, 1, 1), (4, 4, 4)]
       B1.print(); //prints [(1, 1, 1), (6, 6, 6)]
       B2.print(); //prints [(-2, 1, 1), (4, 4, 4)]
       @endcode
    */
    inline Box extrude(const Point& a_dir, int a_dist = 1) const;
    
    /// Extrude (Int Version)
    /**
       AKA single sided grow. Returns a new Box with the <code> a_dir </code> direction extruded a distance <code> a_dist </code).
       By default, the extrusion is upwards. If <code>a_upper</code> is false, extrudes downwards.
        
       \param a_dir    Axis to extrude along
       \param a_dist   (Optional) Distance to extrude (default: 1)
       \param a_upper  (Optional) Extrude upwards? (default: true)

       Example:
       @code
       //DIM=2
       using namespace Proto;
       Box B0 = Box::Cube(4).shift(Point::Ones());
       Box B1 = B0.extrude(2,3,false);
       Box B2 = B0.extrude(2,3,true);
       B0.print(); //prints [(1, 1, 1), (4, 4, 4)]
       B1.print(); //prints [(1, 1, -2), (4, 4, 4)]
       B2.print(); //prints [(1, 1, 1), (4, 4, 7)]
       @endcode
    */
    inline Box extrude(int a_dir, int a_dist = 1, bool a_upper = true) const;

    /// Grow High Side
    /**
        Alias for <code>Box::extrude(idir, igrow, true)</code>.
    */
    inline Box growHi(int idir, int igrow) const;
  
    /// Grow Low Side
    /**
        Alias for <code>Box::extrude(idir, igrow, false)</code>.
    */
    inline Box growLo(int idir, int igrow) const;

    /// Isotropic Coarsen Operation
    /**
       Returns a new Box coarsened by a factor of <code>a_numpoints</code>. Fails if user tries to 
       coarsen using a non-positive ratio.
       If the limits of <code>*this</code> are not multiples of the coarsening ratio (e.g. if <code>this-> coarsenable() != true</code>) the resulting Box is not guaranteed to be a subset of <code>*this</code>. In this situation, use <code>Box::taperCoarsen</code>

       \param a_ratio     Coarsening ratio

       Example:
       @code
       //DIM=2
       using namespace Proto;
       Box B0 = Box::Cube(4);
       Box B1 = B0.coarsen(2);
       B0.print(); //prints [(0, 0), (3, 3)]
       B1.print(); //prints [(0, 0), (1, 1)]

       Box B2 = Box::Cube(3).shift(Point::Ones(2));
       Box B3 = B2.coarsen(2);
       B2.print(); //prints [(2, 2), (4, 4)]
       B3.print(); //prints [(1, 1), (2, 2)]
       @endcode
    */
    inline Box coarsen(unsigned int a_ratio) const;
    
    /// Anisotropic Coarsen Operation
    /**
       Returns a new Box coarsened in each direction according to <code>a_pt</code>.
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


    /// Coarsenable Query
    /**
       Returns true if *this will coarsen normally by a given ratio. 
    
       \param a_ratio    Coarsening ratio
    */
    inline bool coarsenable(const int& a_ratio) const;

    /// Tapered Coarsen
    /**
       This function is identical to Box::coarsen when Box::coarsenable() is true.
       For non-coarsenable Boxes, the new Box's limits are rounded such that the result
       Is always a subset of *this. Specifically, the lower limits are always rounded UP.

       \param a_ref  Coarsening ratio
    */
    inline Box taperCoarsen(const Point& a_ref) const;

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

    /// Edge
    /**
       Returns the subset on the boundary of *this in a given direction with specified thickness.
       Very handy for adjacency and boundary based computations. 
       The returned <code>Box</code> is always a subset of <code>*this</code>

       \param a_dir  Direction of desired edge with respect to the cell center of *this
       \param a_dist Thickness of the output
 
       Examples:
       @code
       //DIM = 2
       using namespace Proto;
       auto B0 = Box::Cube(4).shift(Point::Ones()); //[(1,1), (4,4)]
       auto B1 = B0.edge(Point::Basis(0));          //[(4,1), (4,4)]
       auto B2 = B0.edge(Point::Ones(), 2);         //[(3,3), (4,4)]
       auto B3 = B0.edge(Point::Basis(1,-1), 2);    //[(1,1), (4,2)]
       @endcode
    */
    inline Box edge(const Point& a_dir, int a_dist = 1) const;
    
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
    
    /// Adjacent Cells
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


    /// Lower Adjacent Cells
    /**
        Alias for <code>adjacent(Point::Basis(a_idir), -1)</code>
    */ 
    inline Box adjCellLo(int a_idir, int a_length) const;

    /// Upper Adjacent Cells
    /**
        Alias for <code>adjacent(Point::Basis(a_idir), +1)</code>
    */ 
    inline Box adjCellHi(int a_idir, int a_length) const;

    ///
    inline Box adjCellSide(int a_idir, int a_length, Side::LoHiSide a_sd) const
    {
      Box retval = adjCellHi(a_idir, a_length);
      if(a_sd == Side::Lo)
      {
        retval = adjCellLo(a_idir, a_length);
      }
      return retval;
    }
      
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
    
    /// Print
    /**
       Prints *this using the format <code>[low, high]</code>
    */ 
    inline void print() const; 
    ///@}

    inline size_t linearSize() const
    {
      size_t retval = 0;
      retval += m_low.linearSize();
      retval += m_high.linearSize();
      retval += sizeof(int);
      return retval;
    }

    inline void linearOut(char* a_buf) const
    {
      size_t charsiz = 0;
      char* charbuf = a_buf;

      m_low.linearOut(charbuf);
      charbuf += m_low.linearSize();
      charsiz += m_low.linearSize();

      m_high.linearOut(charbuf);
      charbuf += m_high.linearSize();
      charsiz += m_high.linearSize();

      int* intbuf = (int* ) charbuf;
      *intbuf = m_size;
      charsiz += sizeof(int);
      charbuf += sizeof(int);
    }

    inline void linearIn(const char* a_buf) 
    {
      size_t charsiz = 0;
      const char* charbuf = a_buf;

      m_low.linearIn(charbuf);
      charbuf += m_low.linearSize();
      charsiz += m_low.linearSize();

      m_high.linearIn(charbuf);
      charbuf += m_high.linearSize();
      charsiz += m_high.linearSize();

      int* intbuf = (int* ) charbuf;
      m_size = *intbuf;
      charsiz += sizeof(int);
      charbuf += sizeof(int);
    }
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
