#ifndef PROTO_ARRAY_H
#define PROTO_ARRAY_H
#include <array>
#include "Proto_MemType.H"
namespace Proto {
template<typename T, size_t N, MemType Mem=MEMTYPE_DEFAULT>
class Array {
public:
    ACCEL_DECORATION
    Array();
    ACCEL_DECORATION
    Array(const Array<T,N,Mem> &arr);
    ACCEL_DECORATION
    Array& operator=(const Array<T,N,Mem> &arr);
    Array& operator=(const std::array<T,N> &arr);
    Array(const std::array<T,N> &arr);
    ACCEL_DECORATION
    ~Array();
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
private:
    T* m_data;
};
#include "implem/Proto_ArrayImplem.H"
}
#endif
