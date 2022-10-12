#ifndef PROTO_ARRAY_H
#define PROTO_ARRAY_H
#include <array>
#include "Proto_MemType.H"
namespace Proto {
template<typename T, size_t N, MemType Mem=MEMTYPE_DEFAULT>
class Array {
public:
    ACCEL_DECORATION
    inline Array();

    template<MemType Ory=MEMTYPE_DEFAULT>
    ACCEL_DECORATION
    inline Array(const Array<T,N,Ory> &arr);

    template<MemType Ory=MEMTYPE_DEFAULT>
    ACCEL_DECORATION
    inline Array& operator=(const Array<T,N,Ory> &arr);

    template<MemType Ory=MEMTYPE_DEFAULT>
    ACCEL_DECORATION
    inline Array(const Array<T,N,Ory> &&arr);

    template<MemType Ory=MEMTYPE_DEFAULT>
    ACCEL_DECORATION
    inline Array& operator=(const Array<T,N,Ory> &&arr);

    inline Array& operator=(const std::array<T,N> &arr);

    inline Array(const Point &a_pt);

    inline Array(const std::array<T,N> &arr);

    ACCEL_DECORATION
    inline ~Array();

    ACCEL_DECORATION
    inline void reset();

    ACCEL_DECORATION
    inline void fill(T value);
    
    ACCEL_DECORATION
    inline T* data();
    
    ACCEL_DECORATION
    inline T* data() const;
    
    ACCEL_DECORATION
    inline T& operator[](size_t i);
    
    ACCEL_DECORATION
    inline T& operator[](size_t i) const;

    inline Array<T,N,Mem> operator+(const Array<T,N,Mem>& a_rhs) const; 
    inline Array<T,N,Mem> operator+(T a_rhs) const; 
    inline void operator+=(const Array<T,N,Mem>& a_rhs); 
    inline void operator+=(T a_rhs); 
    
    inline Array<T,N,Mem> operator-(const Array<T,N,Mem>& a_rhs) const; 
    inline Array<T,N,Mem> operator-(T a_rhs) const; 
    inline void operator-=(const Array<T,N,Mem>& a_rhs); 
    inline void operator-=(T a_rhs); 
    
    inline Array<T,N,Mem> operator*(const Array<T,N,Mem>& a_rhs) const; 
    inline Array<T,N,Mem> operator*(T a_rhs) const; 
    inline void operator*=(const Array<T,N,Mem>& a_rhs); 
    inline void operator*=(T a_rhs); 

    // check for division by 0
    inline Array<T,N,Mem> operator/(const Array<T,N,Mem>& a_rhs) const; 
    inline Array<T,N,Mem> operator/(T a_rhs) const; 
    inline void operator/=(const Array<T,N,Mem>& a_rhs); 
    inline void operator/=(T a_rhs); 

    inline T max() const;
    inline T min() const;
    inline T absMax() const;
    inline T sum() const;
    inline T norm() const; //square root of sum of squares

private:
    T* m_data;
};

// Premultiplication by a scalar
template<typename T, unsigned int N, MemType MEM>
inline Array<T,N,MEM> operator*(T a_scale, Array<T,N,MEM> a_array);

// Unary negation
template<typename T, unsigned int N, MemType MEM>
inline Array<T,N,MEM> operator-(Array<T,N,MEM> a_array);

#include "implem/Proto_ArrayImplem.H"
}
#endif
