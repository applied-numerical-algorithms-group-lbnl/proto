#ifndef PROTO_ARRAY_H
#define PROTO_ARRAY_H
#include <array>
namespace Proto {
template<typename T, size_t N>
class Array {
public:
    ACCEL_DECORATION
    inline Array();

    ACCEL_DECORATION
    inline Array(const Array<T,N> &arr);

    ACCEL_DECORATION
    inline Array& operator=(const Array<T,N> &arr);

    ACCEL_DECORATION
    inline Array(const Array<T,N> &&arr);

    ACCEL_DECORATION
    inline Array& operator=(const Array<T,N> &&arr);

    inline Array& operator=(const std::array<T,N> &arr);

    inline Array(const Point &a_pt);

    inline Array(const std::array<T,N> &arr);

    ACCEL_DECORATION
    inline ~Array() = default;

    ACCEL_DECORATION
    inline void reset();

    ACCEL_DECORATION
    inline bool empty();

    ACCEL_DECORATION
    inline void fill(T value);
    
    ACCEL_DECORATION
    inline const T* data();
    
    ACCEL_DECORATION
    inline const T* data() const;
    
    ACCEL_DECORATION
    inline T& operator[](size_t i);
    
    ACCEL_DECORATION
    inline T& operator[](size_t i) const;

    ACCEL_DECORATION
    inline Array<T,N> operator+(const Array<T,N>& a_rhs) const; 
    ACCEL_DECORATION
    inline Array<T,N> operator+(T val) const; 
    ACCEL_DECORATION
    inline void operator+=(const Array<T,N>& a_rhs); 
    ACCEL_DECORATION
    inline void operator+=(T val); 
    
    ACCEL_DECORATION
    inline Array<T,N> operator-(const Array<T,N>& a_rhs) const; 
    ACCEL_DECORATION
    inline Array<T,N> operator-(T val) const; 
    ACCEL_DECORATION
    inline void operator-=(const Array<T,N>& a_rhs); 
    ACCEL_DECORATION
    inline void operator-=(T val); 
    
    ACCEL_DECORATION
    inline Array<T,N> operator*(const Array<T,N>& a_rhs) const; 
    ACCEL_DECORATION
    inline Array<T,N> operator*(T val) const; 
    ACCEL_DECORATION
    inline void operator*=(const Array<T,N>& a_rhs); 
    ACCEL_DECORATION
    inline void operator*=(T val); 

    // check for division by 0
    ACCEL_DECORATION
    inline Array<T,N> operator/(const Array<T,N>& a_rhs) const; 
    ACCEL_DECORATION
    inline Array<T,N> operator/(T val) const; 
    ACCEL_DECORATION
    inline void operator/=(const Array<T,N>& a_rhs); 
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
    inline T norm() const; //square root of sum of squares

    ACCEL_DECORATION
    inline T* begin() const;
    ACCEL_DECORATION
    inline const T* cbegin() const;
    ACCEL_DECORATION
    inline T* end() const;
    ACCEL_DECORATION
    inline const T* cend() const;

private:
    T m_data[N];
};

// Premultiplication by a scalar
template<typename T, unsigned int N>
ACCEL_DECORATION
inline Array<T,N> operator*(T a_scale, Array<T,N> a_array);

// Unary negation
template<typename T, unsigned int N>
ACCEL_DECORATION
inline Array<T,N> operator-(Array<T,N> a_array);

#include "implem/Proto_ArrayImplem.H"
}
#endif
