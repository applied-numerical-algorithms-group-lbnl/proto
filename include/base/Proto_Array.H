/* This class is only necessary because HIP doesn't support STL containers, 
   so std::array can't have any of it's member functions/operators accessed
   on the device.
*/
#pragma once
#ifndef PROTO_ARRAY_H
#define PROTO_ARRAY_H

#include <array>
#include <iostream>
#include <cmath>
#include <limits>
#include "Proto_accel.H"
#include "Proto_Point.H"

namespace Proto {

template<typename T, size_t N>
class ArrayIterator;
template<typename T, size_t N>
class ConstArrayIterator;

template<typename T, size_t N>
class Array {

public:

    ACCEL_DECORATION
    inline static Array<T,N> Zeros();

    ACCEL_DECORATION
    inline static Array<T,N> Ones(T a_scale = 1.0);

    ACCEL_DECORATION
    inline Array();

    ACCEL_DECORATION
    inline Array(T val);

    ACCEL_DECORATION
    inline Array(std::initializer_list<T> list);

    ACCEL_DECORATION
    inline Array(const Array<T,N> &arr);

    ACCEL_DECORATION
    inline Array<T,N>& operator=(const Array<T,N> &arr);

    ACCEL_DECORATION
    inline Array(Array<T,N> &&arr) = default;

    ACCEL_DECORATION
    inline Array<T,N>& operator=(Array<T,N> &&arr) = default;

    ACCEL_DECORATION
    inline Array<T,N>& operator=(const std::array<T,N> &arr);

    ACCEL_DECORATION
    operator std::array<T,N>() const
    {
        std::array<T,N> arr;
        for (int ii = 0; ii < N; ii++)
        {
            arr[ii] = (*this)[ii];
        }
        return arr; 
    }

    ACCEL_DECORATION
    inline Array(const Point &a_pt);

    ACCEL_DECORATION
    inline Array(const std::array<T,N> &arr);

    ACCEL_DECORATION
    inline ~Array() = default;

    ACCEL_DECORATION
    inline void reset();

    ACCEL_DECORATION
    inline void fill(T value);
    
    ACCEL_DECORATION
    inline T* data();
    
    ACCEL_DECORATION
    inline const T* data() const;
    
    ACCEL_DECORATION
    inline T& operator[](size_t i);
    
    ACCEL_DECORATION
    inline const T& operator[](size_t i) const;

    ACCEL_DECORATION
    inline bool operator==(const Array<T,N> &a_rhs) const;
    ACCEL_DECORATION
    inline bool operator!=(const Array<T,N> &a_rhs) const { return !(*this==a_rhs); }

    template<size_t M>
    ACCEL_DECORATION
    inline Array<T,N> operator+(const Array<T,M>& a_rhs) const; 
    ACCEL_DECORATION
    inline Array<T,N> operator+(T val) const; 
    template<size_t M>
    ACCEL_DECORATION
    inline void operator+=(const Array<T,M>& a_rhs); 
    ACCEL_DECORATION
    inline void operator+=(T val); 
    
    template<size_t M>
    ACCEL_DECORATION
    inline Array<T,N> operator-(const Array<T,M>& a_rhs) const; 
    ACCEL_DECORATION
    inline Array<T,N> operator-(T val) const; 
    template<size_t M>
    ACCEL_DECORATION
    inline void operator-=(const Array<T,M>& a_rhs); 
    ACCEL_DECORATION
    inline void operator-=(T val); 
    
    template<size_t M>
    ACCEL_DECORATION
    inline Array<T,N> operator*(const Array<T,M>& a_rhs) const; 
    ACCEL_DECORATION
    inline Array<T,N> operator*(T val) const; 
    template<size_t M>
    ACCEL_DECORATION
    inline void operator*=(const Array<T,M>& a_rhs); 
    ACCEL_DECORATION
    inline void operator*=(T val); 

    template<size_t M>
    ACCEL_DECORATION
    inline Array<T,N> operator/(const Array<T,M>& a_rhs) const; 
    ACCEL_DECORATION
    inline Array<T,N> operator/(T val) const; 
    template<size_t M>
    ACCEL_DECORATION
    inline void operator/=(const Array<T,M>& a_rhs); 
    ACCEL_DECORATION
    inline void operator/=(T val); 

    ACCEL_DECORATION
    inline T max() const;
    ACCEL_DECORATION
    inline T min() const;
    ACCEL_DECORATION
    inline T absMax() const;
    ACCEL_DECORATION
    inline T sum() const;
    ACCEL_DECORATION
    inline T dot(const Array<T,N> &a_rhs) const;
    ACCEL_DECORATION
    inline T product() const;
    ACCEL_DECORATION
    inline T norm() const; //square root of sum of squares

    ACCEL_DECORATION
    inline ArrayIterator<T,N> begin();
    ACCEL_DECORATION
    inline ConstArrayIterator<T,N> begin() const;
    ACCEL_DECORATION
    inline ArrayIterator<T,N> end();
    ACCEL_DECORATION
    inline ConstArrayIterator<T,N> end() const;

    inline void print() const;

private:
    T m_data[N];
};

template<typename T, size_t N>
class ArrayIterator
{
    public:
    inline ArrayIterator(Array<T, N>& a_data, int a_index)
    {
        m_data = &a_data;
        m_index = a_index;
    }

    inline T& operator*() {return (*m_data)[m_index]; } 
    inline bool operator==(const ArrayIterator<T, N>& a_iter) const
    {
        if (a_iter.m_data != m_data) {return false;}
        if (a_iter.m_index != m_index) {return false;}
        return true;
    }
    inline bool operator !=(const ArrayIterator<T, N>& a_iter) const
    {
        return !(*this == a_iter);
    }
    inline ArrayIterator& operator++() {m_index++; return *this; }
    inline ArrayIterator& operator--() {m_index--; return *this; }
    private:
    Array<T, N>* m_data;
    int m_index;
};

template<typename T, size_t N>
class ConstArrayIterator
{
    public:
    inline ConstArrayIterator(const Array<T, N>& a_data, int a_index)
    {
        m_data = &a_data;
        m_index = a_index;
    }

    inline const T& operator*() {return (*m_data)[m_index]; } 
    inline bool operator==(const ConstArrayIterator<T, N>& a_iter) const
    {
        if (a_iter.m_data != m_data) {return false;}
        if (a_iter.m_index != m_index) {return false;}
        return true;
    }
    inline bool operator !=(const ConstArrayIterator<T, N>& a_iter) const
    {
        return !(*this == a_iter);
    }
    inline ConstArrayIterator& operator++() {m_index++; return *this; }
    inline ConstArrayIterator& operator--() {m_index--; return *this; }
    private:
    const Array<T, N>* m_data;
    int m_index;
};

// Premultiplication by a scalar
template<typename T, size_t N>
ACCEL_DECORATION
inline Array<T,N>& operator*(int a_scale, Array<T,N>& a_array);

template<typename T, size_t N>
ACCEL_DECORATION
inline Array<T,N>& operator*(double a_scale, Array<T,N>& a_array);

// Unary negation
template<typename T, size_t N>
ACCEL_DECORATION
inline Array<T,N>& operator-(Array<T,N>& a_array);

template<typename T, size_t N>
inline std::ostream& operator<<(std::ostream& a_os, const Array<T,N>& a_array);

#include "implem/Proto_ArrayImplem.H"
}

#endif
