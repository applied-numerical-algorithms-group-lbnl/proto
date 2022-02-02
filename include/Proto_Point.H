#pragma once
#ifndef _PROTO_POINT_Ht__
#define _PROTO_POINT_Ht__
#include "Proto_cuda.H"
#include "Proto_gpu.H"
#include "Proto_LoHiSide.H"

#include "Proto_PAssert.H"

#include <iostream>
#include <array>

namespace Proto
{

#ifdef PR_HDF5
    class HDF5Handler;
#endif
    /// Integer Valued Vector
    /**
      An element \f$\mathbf{p}\in\mathbb{Z}^{DIM}\f$
    */
    class Point
    {
#ifdef PR_HDF5
        friend class Proto::HDF5Handler;
#endif
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
            Note: Functionality is largely deprecated by the variadic constructor.

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
        CUDA_DECORATION inline explicit Point(vals... args) { unpack(0, args...); }

        /// Copy Constructor
        CUDA_DECORATION inline Point(const Point& a_pt);

        ///@}
        //////////////////////////////////////////////////////////////////////////////////////////
        /** @name Static Methods */
        ///@{

        /// Get Ones
        /**
          Returns the unit Point <code>scale*(1,1,...,1)</code>

          \param a_scale  (Optional) value to scale the vector by (default: 1)
        */
        CUDA_DECORATION inline static Point Ones(int a_scale=1);

        /// Get Unit Point
        /**
          Returns the Point <code>(1,1,...,1)</code>
        */
        CUDA_DECORATION inline static Point Unit() { return Ones(1); }

        /// Get Zeros
        /**
          Returns the zero Point <code>(0,0,...,0)</code>
        */
        CUDA_DECORATION inline static Point Zeros();

        /// Get Zeros (Pedantic Spelling)
        /**
            Identical to Zeros() but with alternative spelling.
        */
        CUDA_DECORATION inline static Point Zeroes() { return Zeros(); }
  
        /// Get Zeros (Alternate Spelling)
        /**
            Identical to Zeros() but with alternative spelling.
        */
        CUDA_DECORATION inline static Point Zero() { return Zeros(); }
    
        /// Get Basis 
        /**
          Returns an axis-aligned Point \f$\bf{p}: \bf{p}[k] = scale*\delta_{dir}^k\f$.

          \param a_dir    A coordinate axis.
          \param a_scale  (Optional) value to scale the unit vector by (default: 1)
          */
        CUDA_DECORATION inline static Point Basis(int a_dir, int a_scale=1);

        /// Get Basis (Side Argument)
        /**
          Returns an axis-aligned Point \f$\bf{p}: \bf{p}[k] = s*\delta_{dir}^k\f$.
          Where \f$ s = 1 \f$ if <code> a_side = Side::Hi </code> and 
          where \f$ s = -1 \f$ if <code> a_side = Side::Lo </code>. 

          \param a_dir    A spatial axis (0 for x, 1 for y, etc.)
          \param a_side   Low or high side. Low amounts to scaling by -1.
        */
        CUDA_DECORATION inline static Point Basis(int a_dir, const Side::LoHiSide a_side);

        ///@}
        //////////////////////////////////////////////////////////////////////////////////////////
        /** @name Accessors */
        ///@{

        /// Get Component (Const)
        /**
          Returns a const reference to the a_index component of *this

          \param a_index  Index in [0,DIM)
          */
        CUDA_DECORATION inline const int&  operator[](unsigned char a_index) const;

        /// Get Component (Non-Const)
        /**
          Returns a mutable reference to the a_index component of *this

          \param a_index  Index in [0,DIM)
          */
        CUDA_DECORATION inline int&  operator[](unsigned char a_index);
            
        ///@}
        //////////////////////////////////////////////////////////////////////////////////////////
        /** @name Operators */
        ///@{

        /// Asignment Operator
        /**
          Performs copy assignment.

          \param a_rhs    Another Point
          */
        CUDA_DECORATION inline Point& operator=(const Point& a_rhs);

        /// Less-Than Operator
        /** 
          Used to define an ordering for placing Points in maps etc.
          Uses Lexical comparison. Note that Chombo lexLT is also lexical
          but uses the dimensions in the opposite way.

          \param a_rhs    Another Point
          */
        CUDA_DECORATION inline bool operator<(const Point& a_rhs) const;

        /// Equality Operator
        /**
          Equality is checked componentwise.   
         
          \param a_pt     Another Point
          */
        CUDA_DECORATION inline bool operator==(const Point& a_pt) const;

        /// Inequality Operator
        /**
          Equality is checked componentwise.   

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
          Creates a new Point by adding nref to each component of *this.

          \param a_nref   An Integer scalar
          */
        inline Point operator+(int a_nref) const;

        ///Scalar Subtraction Operator
        /**
          Creates a new Point by subtracting nref to each component of *this.

          \param a_nref   An Integer scalar
          */
        inline Point operator-(int a_nref) const;

        /// Scalar Multiplication Operator
        /**
          Creates a new Point by multiplying each component of *this by nref.

          \param a_nref   An Integer scalar
          */
        inline Point operator*(int a_nref) const;

        /// Scalar Division Operator
        /**
          Creates a new Point by dividing each component of *this by nref.
          Quotients are rounded down. Division of any component by 0 yields an error. 

          \param a_nref   An Integer scalar
          */
        inline Point operator/(int a_nref) const;

        /// Scalar Modulus Operator
        /**
          Creates a new Point by taking the modulus of each component of *this by nref.
          Quotients are rounded down. Modulo of any component by 0 yields an error. 

          \param a_nref   An Integer scalar
          */
        inline Point operator%(int a_nref) const;

        /// In Place Componentwise Addition Operator
        /**
            Adds another Point componentwise to *this.

            \param a_pt     Another Point
        */
        CUDA_DECORATION inline void operator+=(const Point& a_pt);

        /// In Place Componentwise Subtraction Operator
        /**
            Subtracts another Point componentwise from *this.

            \param a_pt     Another Point
        */
        CUDA_DECORATION inline void operator-=(const Point& a_pt);

        /// In Place Componentwise Multiplication Operator
        /**
            Multiplies *this by another Point componentwise.

            \param a_pt     Another Point
        */
        CUDA_DECORATION  inline void operator*=(const Point& a_pt);

        /// In Place Componentwise Division Operator
        /**
          Divides *this by another Point componentwise.
          Quotients are rounded down. Division of any component by 0 yields an error. 

          \param a_pt       Another Point
          */
        CUDA_DECORATION inline void operator/=(const Point& a_pt);

        /// In Place Componentwise Modulus Operator
        /**
            Computes the modulus of *this by another Point componentwise.
            Quotients are rounded down. Module of any component by 0 yields an error.

            \param a_pt     Another Point
        */
        CUDA_DECORATION inline void operator%=(const Point& a_pt);

        /// In Place Scalar Addition Operator
        /**
            Adds an integer to every component of *this.

            \param a_n  An integer.
        */
        inline void operator+=(int a_n);

        /// In Place Scalar Subtraction Operator
        /**
            Subtracts an integer from every component of *this.

            \param a_n  An integer.
        */
        inline void operator-=(int a_n);

        /// In Place Scalar Multiplication Operator
        /**
            Multiplies every component of *this by an integer.

            \param a_n  An integer.
        */
        inline void operator*=(int a_n);

        /// In Place Scalar Division Operator
        /**
            Divides every component of *this by an integer.
            Quotients are rounded down. Division by 0 yields an error.

            \param a_n  An integer.
        */
        inline void operator/=(int a_n);

        /// In Place Scalar Modulus Operator
        /**
            Computes the modulus of every component of *this by an integer.
            Modulo of 0 results in an error.

            \param a_n  An integer.
        */
        inline void operator%=(int a_n); 
        
        /// Coarsen Operator
        /**
          Returns a new Point coarsened by a_refRatio.
          This function is identical to any of the division operators.

          \param a_refRatio   A non-zero refinement ratio.
        */
        inline Point coarsen(unsigned int a_refRatio) const;
        
        /// Sum Operator
        /**
          Computes the sum of the entries in *this
          */
        inline int sum() const;

        /// Absolute Value Sum Operator
        /**
          Computes the sum of the absolute values of the entries in *this
          */
        inline int absSum() const;

        /// Isotropic Query
        /**
          Checks if all elements of *this are equal.
          */
        inline bool isotropic() const;

        /// Linear Size
        /**
         *  Returns the size of *this in bytes
        */
        inline size_t linearSize() const
        {
            size_t retval = DIM*sizeof(int);
            return retval;
        }

        /// Linear Out
        /**
         *  Writes from *this into a linear buffer of size <code>this->linearSize()</code>.
        */
        inline void linearOut(char* a_buf) const
        {
            int* intbuf = (int *) a_buf;
            for(int idir = 0; idir < DIM; idir++)
            {
                intbuf[idir] = m_tuple[idir];
            }
        }

        /// Linear In
        /**
         *  Reads to *this from a linear buffer of size <code>this->linearSize()</code>.
        */
        inline void linearIn(const char* a_buf) 
        {
            const int* intbuf = (int *) a_buf;
            for(int idir = 0; idir < DIM; idir++)
            {
                m_tuple[idir] = intbuf[idir];
            }
        }

        ///@}
        //////////////////////////////////////////////////////////////////////////////////////////
        /** @name Utility */
        ///@{

        /// Print Function
        /**
          Output is formatted: <code>(p1,p2,...)</code>
          */
        inline void print() const;

        ///@}


        private:

        int m_tuple[DIM]; ///<Integer coordinates of the Point
        
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
    CUDA_DECORATION
    inline void Point::unpack(int D, int i)
    {
        if (D < DIM)
        {
            m_tuple[D] = i;
        }
    }

    //////////////////////////////////////////////////////////////////////////////////////////
    /** @name External Operators */
    ///@{

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

    ///@}

#include "implem/Proto_PointImplem.H"
} //end namespace Proto
#endif //end include guard
