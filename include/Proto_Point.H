#ifndef _PROTO_POINT_Ht__
#define _PROTO_POINT_Ht__
#include "Proto_cuda.H"

#include "Proto_PAssert.H"

#include <iostream>
#include <array>

namespace Proto 
{
    
/// Integer Valued Vector
/**
   An element of \f$\bf{p}\in\bf{Z}^{DIM}\f$
*/
  class Point
  {
  public:
    //////////////////////////////////////////////////////////////////////////////////////////
    /** @name Constructors */
    ///@{
    /// Default Constructor
    /**
       Returns the zero Point <code>(0,0,...,0)</code>
    */
    CUDA_DECORATION inline Point();

    /// C Array Constructor
    /**
       Allows "literal" construction of points.

       Example:
       @code
       using namespace Proto;
       //DIM=2
       Point p = Point({-5,13}); //creates the point (-5,13)
       @endcode

       \param a_tuple  integer array representing a Point
    */
    CUDA_DECORATION inline Point(const int (&a_tuple)[DIM]);

    /// Variadic Constructor
    /**
       Builds a Point from the first DIM arguments given to this. Any trailing arguments are ignored.
       This constructor is very useful for building dimensionally independent test code.

       Example:
       @code
       Point p(1,2,3,4,5,6);
       if (DIM == 1)
       {
       std::cout << p << std::endl;    //prints (1)
       } else if (DIM == 2) {
       std::cout << p << std::endl;    //prints (1, 2)
       } else if (DIM == 3) {
       std::cout << p << std::endl;    //prints (1, 2, 3)
       } //etc.
       @endcode

       \param args     at least DIM integer arguments representing a Point
    */
    template<typename... vals>
    CUDA_DECORATION inline explicit Point(vals... args)
    {
      //PROTO_ASSERT(sizeof...(args)>=DIM,
      //"Point(vals... args) invalid. %lu arguments provided. At least %i are needed.",
      //sizeof...(args),DIM);
      unpack(0, args...);
    }
    /// Copy Constructor
    CUDA_DECORATION inline Point(const Point& a_pt);
        
    ///@}
    //////////////////////////////////////////////////////////////////////////////////////////
    /** @name Static Methods */
    ///@{

    /// Get Basis 
    /**
       Returns an axis-aligned Point \f$\bf{p}: \bf{p}[k] = scale*\delta_{dir}^k\f$.
            
       \param a_dir    A spatial axis (0 for x, 1 for y, etc.)
       \param a_scale  (Optional) value to scale the unit vector by (default: no scaling)
    */
    inline static Point Basis(int a_dir, int a_scale=1);

    ///
    inline int sum() const 
    {
      int retval = 0;
      for(int idir = 0; idir < DIM; idir++)
      {
        retval += (*this)[idir];
      }
      return retval;
    }
    /// Get Zeros
    /**
       Returns the zero Point <code>(0,0,...,0)<code>
    */
      inline static Point Zeros();


    /// Get Zeros (if you like spelling things correctly)
    /**
       Returns the zero Point <code>(0,0,...,0)<code>
    */
    inline static Point Zeroes()
    {
      return Zeros();
    }
  
    // to make chombo conversion a bit saner
    inline static Point Zero()
    {
      return Zeros();
    }

    // to make chombo conversion a bit saner
    inline static Point Unit()
    {
      return Ones(1);
    }

    ///
    inline Point coarsen(unsigned int a_refRat) const
    {
      Point retval = *this;
      for(int idir = 0; idir < DIM; idir++)
      {
        retval[idir] /= a_refRat; 
      }
      return retval;
    }

    /// Get Ones
    /**
       Returns the unit Point <code>scale*(1,1,...,1)</code>

       \param a_scale  (Optional) value to scale the vector by (default: no scaling)
    */
    inline static Point Ones(int a_scale=1);
    ///@}
    //////////////////////////////////////////////////////////////////////////////////////////
    /** @name Accessors */
    ///@{

    /// Get Component (Const)
    /**
       Returns a const reference to the a_index component of *this

       \param a_index  Index in [0,DIM)
    */
    CUDA_DECORATION inline const int&  operator[](unsigned char a_index) const
    {
      PROTO_ASSERT((a_index < DIM),
                   "Point::operator[](uchar a_index) invalid for a_index = %i.\
            a_index must be in [0,DIM = %i)",
                   a_index, DIM);
      return m_tuple[a_index];
    };

    /// Get Component (Non-Const)
    /**
       Returns a mutable reference to the a_index component of *this

       \param a_index  Index in [0,DIM)
    */
    CUDA_DECORATION inline int&  operator[](unsigned char a_index)
    {
      PROTO_ASSERT(((a_index < DIM)),
                   "Point::operator[](uchar a_index) invalid for a_index = %i.\
            a_index must be in [0,DIM = %i)",
                   a_index, DIM);
      return m_tuple[a_index];
    };
            
    ///@}
        
    //////////////////////////////////////////////////////////////////////////////////////////
    /** @name Operators */
    ///@{

    /// Equals Operator
    /**
       Equality is based on value, not instance.

       \param a_rhs    Another Point
    */
    CUDA_DECORATION inline Point& operator=(const Point& a_rhs);

    /// Less-Than Operator
    /** 
        Used to define an ordering for placing Points in maps etc.
        Do not use to actually compare Points.
        Uses Lexical comparison.
        Note that Chombo lexLT is also lexical but uses the dimensions in
        the opposite way.   Both work for maps and so on.
        \param a_rhs    Another Point
    */
    CUDA_DECORATION inline bool operator<(const Point& a_rhs) const;

    /// Equality Operator
    /**
       Returns true if:
       *this[\f$k\f$] == *a_pt[\f$k\f$] \f$\forall\,k\,\in\,[0,DIM)\f$

       \param a_pt     Another Point
    */
    CUDA_DECORATION inline bool operator==(const Point& a_pt) const;

    /// Inequality Operator
    /**
       Returns true if *this != a_pt.

       \param a_pt     Another Point
    */
    CUDA_DECORATION inline bool operator!=(const Point& a_pt) const;
        
    /// Componentwise Point Addition Operator
    /**
       Creates the new Point by adding <code>*this</code> and <code>rhs</code> componentwise.

       \param a_rhs     Another Point
    */
    CUDA_DECORATION inline Point operator+(const Point& a_rhs) const;

    /// Componentwise Point Subtraction Operator 
    /**
       Creates the new Point by subtracting <code>*this</code> 
       and <code>rhs</code> componentwise.

       \param a_rhs     Another Point
    */
    CUDA_DECORATION inline Point operator-(const Point& a_rhs) const;

    /// Componentwise Point Multiplication Operator
    /**
       Creates the new Point by multiplying <code>*this</code> 
       and <code>rhs</code> componentwise.

       \param a_pt     Another Point
    */
    CUDA_DECORATION inline Point operator*(const Point& a_pt) const;

    /// Componentwise Point Division Operator
    /**
       Creates the new Point by dividing <code>*this</code> 
       by <code>rhs</code> componentwise. Quotients are rounded down.
       Division of any component by 0 yields an error. 

       \param a_pt     Another Point
    */
    CUDA_DECORATION inline Point operator/(const Point& a_pt) const;
        
    /// Componentwise Point Modulus Operator
    /**
       Creates the new Point by taking the modulus of  <code>*this</code> 
       by <code>rhs</code> componentwise. Quotients are rounded down.
       Modulo by 0 yields an error. 

       \param a_pt     Another Point
    */
    CUDA_DECORATION inline Point operator%(const Point& a_pt) const;

    ///Scalar Addition Operator
    /**
       Creates a new point by adding nref to each component of *this.

       \param a_nref   An Integer scalar
    */
    inline Point operator+(int a_nref) const;

    ///Scalar Subtraction Operator
    /**
       Creates a new point by subtracting nref to each component of *this.

       \param a_nref   An Integer scalar
    */
    inline Point operator-(int a_nref) const;

    /// Scalar Multiplication Operator
    /**
       Creates a new point by multiplying each component of *this by nref.

       \param a_nref   An Integer scalar
    */
    inline Point operator*(int a_nref) const;

    /// Scalar Division Operator
    /**
       Creates a new point by dividing each component of *this by nref.
       Quotients are rounded down. Division of any component by 0 yields an error. 

       \param a_nref   An Integer scalar
    */
    inline Point operator/(int a_nref) const;

    /// Scalar Modulus Operator
    /**
       Creates a new point by taking the modulus of each component of *this by nref.
       Quotients are rounded down. Modulo of any component by 0 yields an error. 

       \param a_nref   An Integer scalar
    */
    inline Point operator%(int a_nref) const;
        
    /// In Place Componentwise Addition Operator
    CUDA_DECORATION inline void operator+=(const Point& a_pt);

    /// In Place Componentwise Subtraction Operator
    CUDA_DECORATION inline void operator-=(const Point& a_pt);

    /// In Place Componentwise Multiplication Operator
    CUDA_DECORATION  inline void operator*=(const Point& a_pt);

    /// In Place Componentwise Division Operator
    /**
       Quotients are rounded down. Division of any component by 0 yields an error. 

       \param a_pt   Another Point
    */
    CUDA_DECORATION inline void operator/=(const Point& a_pt);

    /// In Place Componentwise Modulus Operator
    CUDA_DECORATION inline void operator%=(const Point& a_pt);
        
    /// In Place Scalar Addition Operator
    inline void operator+=(int a_n);

    /// In Place Scalar Subtraction Operator
    inline void operator-=(int a_n);

    /// In Place Scalar Multiplication Operator
    inline void operator*=(int a_n);

    /// In Place Scalar Division Operator
    /**
       Quotients are rounded down. Division of any component by 0 yields an error. 
    */
    inline void operator/=(int a_n);

    /// In Place Scalar Modulus Operator
    inline void operator%=(int a_n); 

    /// Computes the sum of the entries in *this
    inline int sum()
    {
        int s = 0;
        for (int ii = 0; ii < DIM; ii++) {
            s += m_tuple[ii];
        }
        return s;
    }

    /// Computes the sum of the absolute values of the entries in *this
    inline int absSum()
    {
        int s = 0;
        for (int ii = 0; ii < DIM; ii++) {
            s += std::abs(m_tuple[ii]);
        }
        return s;
    }
    ///@}
    //////////////////////////////////////////////////////////////////////////////////////////
    /** @name Utility */
    ///@{

    /// Print Function
    /**
      Output is formatted: <code>(p1,p2,...)</code>
     */
    inline void print() const
    {
        using std::cout;
        using std::endl;
        cout << "(";
        for (int dir = 0; dir < DIM ; dir++)
        {
            cout << m_tuple[dir];
            if(dir < DIM-1) cout << ",";
        }
        cout <<   ")"  << endl;
    };

    ///@}

    int m_tuple[DIM]; ///<Integer coordinates of the Point

    private:

    template<typename... vars>
        CUDA_DECORATION inline void unpack(int D, int i, vars... args)
        {
            if (D < DIM)
            {
                m_tuple[D] = i;
                unpack(D+1, args...);
            }
        }
  }; //end class Point

  template<>
      inline void Point::unpack(int D, int i)
      {
          if (D < DIM)
          {
              m_tuple[D] = i;
          }
      }

  /// Stream Operator
  inline std::ostream& operator <<(std::ostream& a_os, const Point& a_pt)
  {
      a_os << "(" ;
      for (int dir = 0; dir < DIM; dir++) 
      {
          a_os << a_pt[dir];
          if(a_pt[dir] < 10)
          {
            a_os << " ";
          }
          if (dir < DIM-1) 
          {
              a_os << ",";
          }
      }
      a_os << ")";
      return a_os;
  }

  /// Premultiplication by scalar
  inline Point operator*(int a_scale, Point a_pt)
  {
      return (a_pt*a_scale);
  }

  /// Unary Negation
  inline Point operator-(Point a_pt)
  {
      return a_pt*(-1);
  }

#include "implem/Proto_PointImplem.H"
} //end namespace Proto
#endif //end include guard
