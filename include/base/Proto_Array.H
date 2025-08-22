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

#include "Proto_Accel.H"
#include "Proto_Point.H"

namespace Proto {

template<typename T, size_t N>
class ArrayIterator;
template<typename T, size_t N>
class ConstArrayIterator;

/// @brief A templated constant size array object similar to std::array, but with the ability to be used inside of device kernels
/// @tparam T a plain-old-data type (int, double, etc)
/// @tparam N static size of the array
template<typename T, size_t N>
class Array {

public:

    /// @brief Shortcut for creating an Array filled with zeros
    /// @return An Array of zeros
    ACCEL_DECORATION
    inline static Array<T,N> Zeros();

    /// @brief Shortcut for creating an isotropic array
    /// @param scale - An optional scale. The output vector will be [scale, scale, scale, ...]. [Default: 1]
    /// @return An isotropic Array multiplied by scale
    ACCEL_DECORATION
    inline static Array<T,N> Ones(T scale = 1.0);

    /// @brief Constructs an Array with uninitialized values
    ACCEL_DECORATION
    inline Array();

    /// @brief Constructs an Array with all values set to initValue
    /// @param initValue - The initial value of all elements of the Array
    ACCEL_DECORATION
    inline Array(T initValue);

    /// @brief Constructs an Array using initializer_list syntax
    /** 
     *  Example:
     *  @code
     *  Array<int, 5> a = {1,2,3,4,5};
     *  @endcode
     */
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

    /// @brief Implicit Cast to std::array
    /// @return 
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

    /// @brief Constructs an int-valued Array from a Point
    ACCEL_DECORATION
    inline Array(const Point& point);

    /// @brief Constructs an Array from a std::array
    ACCEL_DECORATION
    inline Array(const std::array<T,N>& arr);

    ACCEL_DECORATION
    inline ~Array() = default;

    // Deprecated
    // ACCEL_DECORATION
    // inline void reset();

    /// @brief Sets all values of *this to value
    ACCEL_DECORATION
    inline void fill(T value);
    
    /// @brief Returns a raw pointer to the data buffer of *this
    ACCEL_DECORATION
    inline T* data();
    
    /// @brief Returns a const raw pointer to the data buffer of *this
    ACCEL_DECORATION
    inline const T* data() const;
    
    /// @brief Returns a reference to the element at index i
    ACCEL_DECORATION
    inline T& operator[](size_t i);
    
    /// @brief Returns a const reference to the element at index i
    ACCEL_DECORATION
    inline const T& operator[](size_t i) const;

    /// @brief Equality operator. Defined by value.
    ACCEL_DECORATION
    inline bool operator==(const Array<T,N> &rhs) const;

    /// @brief Inequality operator. Defined by value.
    ACCEL_DECORATION
    inline bool operator!=(const Array<T,N> &rhs) const { return !(*this==rhs); }

    /// @brief Vector addition operator. Returns a new Array of the same size as *this.
    /// if M < N, the elements of rhs are added to the first M elements of the output
    /// @tparam M - The number of elements in rhs
    /// @param rhs - Another Array 
    /// @return The elementwise sum of *this with rhs on the first M elements
    template<size_t M>
    ACCEL_DECORATION
    inline Array<T,N> operator+(const Array<T,M>& rhs) const; 

    /// @brief Scalar addition operator. Returns a new Array containing the values of *this + value
    ACCEL_DECORATION
    inline Array<T,N> operator+(T value) const;

    /// @brief Vector in-place addition operator. The first M elements of rhs are added to this.
    /// @tparam M - The number of elements in rhs
    /// @param rhs - Another Array 
    template<size_t M>
    ACCEL_DECORATION
    inline void operator+=(const Array<T,M>& rhs); 

    /// @brief Scalar in-place addition operator. Adds value to each element of *this
    ACCEL_DECORATION
    inline void operator+=(T value); 
    
    /// @brief Vector subtraction operator. Returns a new Array of the same size as *this.
    /// if M < N, the elements of rhs are subtracted from the first M elements of the output
    /// @tparam M - The number of elements in rhs
    /// @param rhs - Another Array 
    /// @return The elementwise difference of *this with rhs on the first M elements
    template<size_t M>
    ACCEL_DECORATION
    inline Array<T,N> operator-(const Array<T,M>& rhs) const; 

    /// @brief Scalar subtraction operator. Subtracts value from each element of *this
    ACCEL_DECORATION
    inline Array<T,N> operator-(T value) const; 

    /// @brief Vector in-place subtraction operator. The first M elements of rhs are subtracted from this.
    /// @tparam M - The number of elements in rhs
    /// @param rhs - Another Array 
    template<size_t M>
    ACCEL_DECORATION
    inline void operator-=(const Array<T,M>& rhs); 

    /// @brief Scalar in-place subtraction operator. Subtracts value from each element of *this
    ACCEL_DECORATION
    inline void operator-=(T value); 
    
    /// @brief Vector product operator. Returns a new Array of the same size as *this.
    /// if M < N, the elements of rhs are multiplied with the first M elements of the output
    /// @tparam M - The number of elements in rhs
    /// @param rhs - Another Array 
    /// @return The elementwise product of *this with rhs on the first M elements
    template<size_t M>
    ACCEL_DECORATION
    inline Array<T,N> operator*(const Array<T,M>& rhs) const; 

    /// @brief Scalar multiplication operator. Multiplies each element of *this by value
    ACCEL_DECORATION
    inline Array<T,N> operator*(T value) const; 

    /// @brief Vector in-place multiplication operator. The first M elements of rhs are multplied into this.
    /// @tparam M - The number of elements in rhs
    /// @param rhs - Another Array 
    template<size_t M>
    ACCEL_DECORATION
    inline void operator*=(const Array<T,M>& rhs); 

    /// @brief Scalar in-place multiplication operator. Multiplies each value of *this by value.
    ACCEL_DECORATION
    inline void operator*=(T value); 

    /// @brief Vector quotient operator. Returns a new Array of the same size as *this.
    /// if M < N, the elements of rhs are divided by the first M elements of the output.
    /// Fails by assertion if any element of rhs contains zero
    /// @tparam M - The number of elements in rhs
    /// @param rhs - Another Array 
    /// @return The elementwise quotient of *this with rhs on the first M elements
    template<size_t M>
    ACCEL_DECORATION
    inline Array<T,N> operator/(const Array<T,M>& nonZeroRhs) const; 

    /// @brief Scalar quotient operator. Divides each element of *this by value.
    /// Fails by assertion if value is zero.
    ACCEL_DECORATION
    inline Array<T,N> operator/(T nonZeroValue) const; 

    /// @brief Vector in-place quotient operator. The first M elements this are divided by rhs elementwise
    /// Fails by assertion if any element of rhs is zero.
    /// @tparam M - The number of elements in rhs
    /// @param rhs - Another Array 
    template<size_t M>
    ACCEL_DECORATION
    inline void operator/=(const Array<T,M>& nonZeroRhs); 

    /// @brief Scalar in-place quotient operator. Divides each value of *this by value.
    /// Fails by assertion of value is zero.
    ACCEL_DECORATION
    inline void operator/=(T nonZeroValue); 

    /// @brief Returns the maximum value of this
    ACCEL_DECORATION
    inline T max() const;

    /// @brief Returns the minimum value of this
    ACCEL_DECORATION
    inline T min() const;

    /// @brief Returns the maximum value of the absolute value of this
    ACCEL_DECORATION
    inline T absMax() const;

    /// @brief Returns the sum of the elements of this
    ACCEL_DECORATION
    inline T sum() const;

    /// @brief Returns the dot product of this with rhs.
    /// The dot product is defined as the sum of the elementwise products of the elements of the two Arrays.
    ACCEL_DECORATION
    inline T dot(const Array<T,N> &rhs) const;

    /// @brief Returns the product of the elements of this.
    ACCEL_DECORATION
    inline T product() const;

    /// @brief Returns the L2 norm of the elements of this.
    /// The L2 norm is defined as the square root of the sum of the squared elements of this
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

    /// @brief Returns a string representation of this
    inline std::string str() const;

    /// @brief Print a representation of this to a std::ostream [Default: std::cout]
    inline void print(std::ostream& stream = std::cout) const;

private:
    T m_data[N];
};

template<typename T, size_t N>
class ArrayIterator
{
    public:
    inline ArrayIterator(Array<T, N>& arr, int index)
    {
        m_data = &arr;
        m_index = index;
    }

    inline T& operator*() {return (*m_data)[m_index]; } 
    inline bool operator==(const ArrayIterator<T, N>& iter) const
    {
        if (iter.m_data != m_data) {return false;}
        if (iter.m_index != m_index) {return false;}
        return true;
    }
    inline bool operator !=(const ArrayIterator<T, N>& iter) const
    {
        return !(*this == iter);
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
    inline ConstArrayIterator(const Array<T, N>& arr, int index)
    {
        m_data = &arr;
        m_index = index;
    }

    inline const T& operator*() {return (*m_data)[m_index]; } 
    inline bool operator==(const ConstArrayIterator<T, N>& iter) const
    {
        if (iter.m_data != m_data) {return false;}
        if (iter.m_index != m_index) {return false;}
        return true;
    }
    inline bool operator !=(const ConstArrayIterator<T, N>& iter) const
    {
        return !(*this == iter);
    }
    inline ConstArrayIterator& operator++() {m_index++; return *this; }
    inline ConstArrayIterator& operator--() {m_index--; return *this; }
    private:
    const Array<T, N>* m_data;
    int m_index;
};

/// @brief Premultiplication by a scalar int
template<typename T, size_t N>
ACCEL_DECORATION
inline Array<T,N>& operator*(int scale, Array<T,N>& arr);

/// @brief Premultiplication by a scalar double
template<typename T, size_t N>
ACCEL_DECORATION
inline Array<T,N>& operator*(double scale, Array<T,N>& arr);

/// @brief Unary negation
template<typename T, size_t N>
ACCEL_DECORATION
inline Array<T,N>& operator-(Array<T,N>& arr);

/// @brief Ostream operator
template<typename T, size_t N>
inline std::ostream& operator<<(std::ostream& stream, const Array<T,N>& arr);

}

#endif
