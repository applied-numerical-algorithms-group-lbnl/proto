#pragma once
#ifndef _PROTO_MATRIX_H_
#define _PROTO_MATRIX_H_

#include "Proto.H"
//TODO: hook this up into the build system
#include "blis.h"
#include "Proto_Lapack.H"

namespace Proto {

template<typename T> class MatrixValue;

template<typename T>
class Matrix
{
    public:
   
    inline static Matrix<T> Identity();

    /// Default Constructor
    inline Matrix(){};
    
    /// Size Constructor
    inline Matrix(unsigned int a_numRows, unsigned int a_numCols);
    
    /// Construct from Existing Buffer
    /** Elements are assumed to be stored in column-major order.
     * 
     * \param a_buffer  elements of the Matrix in column-major order
     * \param a_numRows number of rows
     * \param a_numCols number of columns
    */
    inline Matrix(std::shared_ptr<T>& a_buffer, unsigned int a_numRows, unsigned int a_numCols);

    /// Move Constructor
    inline Matrix(const Matrix<T>&& a_matrix);

    /// Initialization Constructor
    inline Matrix(std::initializer_list<std::initializer_list<T>> a_args);

    /// Destructor
    inline ~Matrix();

    /// Allocate
    inline void define(unsigned int a_numRows, unsigned int a_numCols);

    /// Copy Assignment (Disabled)
    inline Matrix<T>& operator=(const Matrix<T>& a_matrix) = delete;

    /// Move Assignment
    inline Matrix<T>& operator=(const Matrix<T>&& a_matrix);

    /// Equality
    inline bool operator==(const Matrix<T>& a_matrix) const;
    
    /// Inequality
    inline bool operator!=(const Matrix<T>& a_matrix) const;

    /// Get Number of Rows
    inline unsigned int numRows() const;

    /// Get Number of Columns
    inline unsigned int numCols() const;

    /// Get Number of Rows (Shorthand)
    inline unsigned int M() const {return numRows();}
    
    /// Get Number of Columns (Shorthand)
    inline unsigned int N() const {return numCols();}

    /// Row Alias
    /** Create a new Matrix which is an alias to a specified row of *this.
     *
     * \param a_rowIndex    The index of the desired row
     */
    inline Matrix<T> row(unsigned int a_rowIndex);
    
    /// Column Alias
    /** Create a new Matrix which is an alias to a specified row of *this.
     *
     * \param a_rowIndex    The index of the desired row
     */
    inline Matrix<T> col(unsigned int a_colIndex);
    
    /// General Alias
    /** Create a new Matrix which aliases a subset of *this. The alias is the
     * numRows X numCols Matrix with origin at (rowOffset, colOffset).
     *
     * \param a_numRows     Number of rows in the alias
     * \param a_numCols     Number of columns in the alias
     * \param a_rowOffset   Index of the first row of the alias
     * \param a_colOffset   Index of the first column of the alias
     */
    inline Matrix<T> slice(
            unsigned int a_numRows, unsigned int a_numCols,
            unsigned int a_rowOffset, unsigned int a_colOffset);
    inline const Matrix<T> slice(
            unsigned int a_numRows, unsigned int a_numCols,
            unsigned int a_rowOffset, unsigned int a_colOffset) const;

    /// Get Value
    /** Retrieve the value at indices (i, j). At the moment, this operation
     *  cannot be done by reference.
     *
     *  \param a_i  Row index
     *  \param a_j  Column index
     */
    inline T get(unsigned int a_i, unsigned int a_j) const;
    
    /// Get Value (Linear Index)
    /** Retrieve the value at linear index k. The linear index is based on 
     *  column-major ordering: k = colIndex*numRows + rowIndex. At the moment, this operation
     *  cannot be done by reference.
     *
     *  \param a_k  Linear index
     */
    inline T get(unsigned int a_k) const;
   
    inline MatrixValue<T> operator()(unsigned int a_i, unsigned int a_j) { return MatrixValue<T>(*this, a_i, a_j); }
    inline MatrixValue<T> operator()(unsigned int a_k) { return MatrixValue<T>(*this, a_k % M(), a_k / M()); }
    inline const MatrixValue<T> operator() (unsigned int a_i, unsigned int a_j) const { return MatrixValue<T>(*this, a_i, a_j); }
    inline const MatrixValue<T> operator() (unsigned int a_k) const { return MatrixValue<T>(*this, a_k % M(), a_k / M()); }

    /// Set Value
    /** Set the value at indices (i, j)
     *
     *  \param a_i      Row index
     *  \param a_j      Column index
     *  \param a_val    Assigned value
     */
    inline void set(unsigned int a_i, unsigned int a_j, const T& a_val);
    
    /// Set Value (Linear index)
    /** Set the value at linear index k. The linear index is based on 
     *  column-major ordering: k = colIndex*numRows + rowIndex.
     *
     *  \param a_k      Linear index
     *  \param a_val    Assigned value
     */
    inline void set(unsigned int a_k, const T& a_val);
    
    /// Set All Values
    inline void set(const T& a_val);
    inline void setVal(const T& a_val) { set(a_val); }
  
    /// Copy 
    inline void copyTo(Matrix<T>& a_rhs) const;

    /// Add Two Matrices
    /** The left and right hand sides must have the same dimensions
     */
    inline Matrix<T> operator+(const Matrix<T>& a_rhs) const;
    
    /// Add Two Matrices (In Place)
    /** The left and right hand sides must have the same dimensions
     */
    inline void operator+=(const Matrix<T>& a_rhs);
    
    /// Add Constant
    inline Matrix<T> operator+(T a_scale) const;
    
    /// Add Constant (In Place)
    inline void operator+=(T a_scale);
    
    /// Subtract Two Matrices
    /** The left and right hand sides must have the same dimensions
     */
    inline Matrix<T> operator-(const Matrix<T>& a_rhs) const;
    
    /// Add Two Matrices (In Place)
    /** The left and right hand sides must have the same dimensions
     */
    inline void operator-=(const Matrix<T>& a_rhs);
   
    /// Subtract Two Matrices
    inline Matrix<T> operator-(T a_scale) const;
    
    /// Subtract Two Matrices (In Place)
    inline void operator-=(T a_scale);

    /// Multiply Two Matrices
    /** Computes the product (*this)*a_rhs. this->numCols() must be equal
     *  to a_rhs.numRows(). The result is a this->numRows() X a_rhs.numCols() matrix.
     *  This function is the equivalent to (and indeed is a wrapper for) LAPACK's dgemm.
     */
    inline Matrix<T> operator*(const Matrix<T>& a_rhs) const;
    
    /// Scale a Matrix
    inline Matrix<T> operator*(T a_scale) const;
    
    /// Scale a Matrix (In Place)
    inline void operator*=(T a_scale);
    
    /// Tensor Contraction
    /** Returns the sum of element-wise products of two matrices. If both objects are vectors of
     *  the same length, orientation is ignored and the dot product is returned.
     */
    inline T dot(const Matrix<T>& a_rhs) const;

    inline T* data();
    inline const T* data() const;

    /// Transpose
    /** Returns an read-only alias to *this which has been transposed. If a mutable transpose
     *  is needed, use copyTo.
     */
    inline const Matrix<T> transpose() const;

    inline double conditionNumber() const;

    inline void checkConditionNumber() const;

    /// Inverse
    /**
     */
    inline Matrix<T> inverse() const;

    /// Size
    inline unsigned int size() const;

    /// Size (in bytes)
    inline size_t linearSize() const;

    /// Print
    inline void print(std::string a_format = "") const;
    
    /// Vector Query
    /** Returns true if either dimension of *this is equal to 1 */
    inline bool isVector() const;

    /// Scalar Query
    /** Returns true if both dimensions of *this are equal to 1 */
    inline bool isScalar() const;

    private:

    inline num_t blisDatatype() const {}

    bool  m_alias = false;
    obj_t m_object;
};

template<>
inline num_t Matrix<double>::blisDatatype() const {return BLIS_DOUBLE;}

template<typename T>
class MatrixValue
{
    public:

    inline MatrixValue(Matrix<T>& a_matrix, unsigned int a_i, unsigned int a_j)
    : m_matrix(&a_matrix), m_i(a_i), m_j(a_j) {}

    inline T& operator=(T a_value) { m_matrix->set(m_i, m_j, a_value); return a_value;}
    
    inline operator double() const { return (double)m_matrix->get(m_i, m_j); }

    template<typename U>
    inline bool operator==(const U& a_rhs) const { return ((double)(*this) == (double)a_rhs); }
    template<typename U>
    inline bool operator!=(const U& a_rhs) const { return !((*this) == a_rhs); }
    private:
        Matrix<T>* m_matrix;
        unsigned int m_i;
        unsigned int m_j;
};

template<typename T>
inline void solveLS(Matrix<T>& a_A, Matrix<T>& a_b);

template<typename T>
inline void solveRRLS(Matrix<T>& a_A, Matrix<T>& a_b);

#include "implem/Proto_MatrixImplem.H"
} //end namespace Proto
#endif //end include guard
