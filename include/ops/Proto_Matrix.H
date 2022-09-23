#pragma once
#ifndef _PROTO_MATRIX_H_
#define _PROTO_MATRIX_H_

//TODO: hook this up into the build system
#include "../../blis/include/sandybridge/blis.h"

namespace Proto {

template<typename T>
class Matrix
{
    public:
   
    inline static Matrix<T> Identity();

    inline Matrix(){};
    inline Matrix(unsigned int a_numRows, unsigned int a_numCols);
    inline Matrix(std::shared_ptr<T> a_buffer, unsigned int a_numRows, unsigned int a_numCols);
    inline Matrix(const Matrix<T>&& a_matrix);
    inline ~Matrix();
    inline Matrix<T>& operator=(const Matrix<T>& a_matrix) = delete;
    inline Matrix<T>& operator=(const Matrix<T>&& a_matrix);

    inline unsigned int numRows() const;
    inline unsigned int numCols() const;
    inline Matrix<T> row(unsigned int a_rowIndex);
    inline Matrix<T> col(unsigned int a_colIndex);
    inline Matrix<T> slice(
            unsigned int a_numRows, unsigned int a_numCols,
            unsigned int a_rowOffset, unsigned int a_colOffset);
    inline Matrix<T> operator*(const Matrix<T>& a_rhs) const;
    inline void operator*=(const Matrix<T>& a_rhs);
    inline Matrix<T> operator*(T a_scale) const;
    inline void operator*=(T a_scale);
    
    inline Matrix<T> operator+(const Matrix<T>& a_rhs) const;
    inline void operator+=(const Matrix<T>& a_rhs);
    inline Matrix<T> operator+(T a_scale) const;
    inline void operator+=(T a_scale);
    
    inline T dot(const Matrix<T>& a_rhs) const;
    
    inline void print() const;
    inline bool isVector() const;

    private:

    inline num_t blisDatatype(){}

    bool  m_alias = false;
    obj_t m_object;
};

template<>
inline num_t Matrix<double>::blisDatatype() {return BLIS_DOUBLE;}

#include "implem/Proto_MatrixImplem.H"
} //end namespace Proto
#endif //end include guard
