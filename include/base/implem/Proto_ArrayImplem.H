#include "../Proto_Array.H"

namespace Proto {

template<typename T, size_t N>
ACCEL_DECORATION
Array<T,N> Array<T,N>::Zeros()
{
    return Array(0);
}

template<typename T, size_t N>
ACCEL_DECORATION
Array<T,N> Array<T,N>::Ones(T a_scale)
{
    return Array(a_scale);
}

template<typename T, size_t N>
ACCEL_DECORATION
Array<T,N>::Array() {}

template<typename T, size_t N>
ACCEL_DECORATION
Array<T,N>::Array(T val) {
    for (size_t ii = 0; ii < N; ii++) { m_data[ii] = val; } 
}

template<typename T, size_t N>
ACCEL_DECORATION
Array<T,N>::Array(std::initializer_list<T> list) {
    size_t ii = 0;
    for (auto& item : list)
    {
        m_data[ii] = item;
        ii++;
        if (ii >= N) {break;}
    }
} 

template<typename T, size_t N>
ACCEL_DECORATION
Array<T,N>::Array(const Array<T,N> &arr) {
    for (size_t ii = 0; ii < N; ii++) { m_data[ii] = arr[ii]; }
}

template<typename T, size_t N>
ACCEL_DECORATION
Array<T,N>& Array<T,N>::operator=(const Array<T,N> &arr) {
    for (size_t ii = 0; ii < N; ii++) { m_data[ii] = arr[ii]; }
    return *this;
}

template<typename T, size_t N>
ACCEL_DECORATION
Array<T,N>& Array<T,N>::operator=(const std::array<T,N> &arr) {
    for (size_t ii = 0; ii < N; ii++) { m_data[ii] = arr[ii]; }
    return *this;
}

template<typename T, size_t N>
ACCEL_DECORATION
Array<T,N>::Array(const std::array<T,N> &arr) {
    for (size_t ii = 0; ii < N; ii++) { m_data[ii] = arr[ii]; }
}

template<typename T, size_t N>
ACCEL_DECORATION
Array<T,N>::Array(const Point &a_pt) {
    for (int i=0; i<std::min(DIM,(int)N); i++) { m_data[i] = (T)a_pt[i]; }
}

template<typename T, size_t N>
ACCEL_DECORATION
bool Array<T,N>::operator==(const Array<T,N> &a_rhs) const {
    for (size_t i=0; i<N; i++) {
        if (m_data[i] != a_rhs[i]) { return false; }
    }
    return true;
}

template<typename T, size_t N>
template<size_t M>
ACCEL_DECORATION
Array<T,N> Array<T,N>::operator+(const Array<T,M> &a_rhs) const {
    Array<T,N> ret = *this;
    for (size_t i=0; i<std::min(N,M); i++) { ret[i] += a_rhs[i]; }
    return ret;
}

template<typename T, size_t N>
ACCEL_DECORATION
Array<T,N> Array<T,N>::operator+(T val) const {
    Array<T,N> ret = *this;
    for (size_t ii = 0; ii < N; ii++) { ret[ii] += val; }
    return ret;
}

template<typename T, size_t N>
template<size_t M>
ACCEL_DECORATION
void Array<T,N>::operator+=(const Array<T,M> &a_rhs) {
    for (size_t i=0; i<std::min(N,M); i++) { m_data[i] += a_rhs[i]; }
}

template<typename T, size_t N>
ACCEL_DECORATION
void Array<T,N>::operator+=(T val) {
    for (size_t ii = 0; ii < N; ii++) {m_data[ii] += val; }
}

template<typename T, size_t N>
template<size_t M>
ACCEL_DECORATION
Array<T,N> Array<T,N>::operator-(const Array<T,M> &a_rhs) const {
    Array<T,N> ret = *this;
    for (size_t i=0; i<std::min(N,M); i++) { ret[i] -= a_rhs[i]; }
    return ret;
}

template<typename T, size_t N>
ACCEL_DECORATION
Array<T,N> Array<T,N>::operator-(T val) const {
    Array<T,N> ret = *this;
    for (size_t ii = 0; ii < N; ii++) { ret[ii] -= val; }
    return ret;
}

template<typename T, size_t N>
template<size_t M>
ACCEL_DECORATION
void Array<T,N>::operator-=(const Array<T,M> &a_rhs) {
    for (size_t i=0; i<std::min(N,M); i++) { m_data[i] -= a_rhs[i]; }
}

template<typename T, size_t N>
ACCEL_DECORATION
void Array<T,N>::operator-=(T val) {
    for (size_t ii = 0; ii < N; ii++) { m_data[ii] -= val; }
}

template<typename T, size_t N>
template<size_t M>
ACCEL_DECORATION
Array<T,N> Array<T,N>::operator*(const Array<T,M> &a_rhs) const {
    Array<T,N> ret = *this;
    for (size_t ii = 0; ii < N; ii++) { ret[ii] *= a_rhs[ii]; }
    return ret;
}

template<typename T, size_t N>
ACCEL_DECORATION
Array<T,N> Array<T,N>::operator*(T val) const {
    Array<T,N> ret = *this;
    for (size_t ii = 0; ii < N; ii++) { ret[ii] *= val; }
    return ret;
}

template<typename T, size_t N>
template<size_t M>
ACCEL_DECORATION
void Array<T,N>::operator*=(const Array<T,M> &a_rhs) {
    for (size_t i=0; i<std::min(M,N); i++) { m_data[i] *= a_rhs[i]; }
}

template<typename T, size_t N>
ACCEL_DECORATION
void Array<T,N>::operator*=(T val) {
    for (size_t ii = 0; ii < N; ii++) {m_data[ii] *= val; }
}

// check for division by 0
template<typename T, size_t N>
template<size_t M>
ACCEL_DECORATION
Array<T,N> Array<T,N>::operator/(const Array<T,M> &a_rhs) const {
    Array<T,N> ret = *this;
    for (size_t i=0; i<std::min(M,N); i++) {
        PR_assert(a_rhs[i]!=0)
        ret[i] /= a_rhs[i];
    }
    return ret;
}

template<typename T, size_t N>
ACCEL_DECORATION
Array<T,N> Array<T,N>::operator/(T val) const {
    PR_assert(val!=0)
    Array<T,N> ret = *this;
    for (size_t ii = 0; ii < N; ii++) {ret[ii] /= val; }
    return ret;
}

template<typename T, size_t N>
template<size_t M>
ACCEL_DECORATION
void Array<T,N>::operator/=(const Array<T,M> &a_rhs) {
    for (size_t i=0; i<std::min(M,N); i++) {
        PR_assert(a_rhs[i]!=0)
        m_data[i] /= a_rhs[i];
    }
}

template<typename T, size_t N>
ACCEL_DECORATION
void Array<T,N>::operator/=(T val) {
    PR_assert(val!=0)
    for (size_t ii = 0; ii < N; ii++) {m_data[ii] /= val; }
}

template<typename T, size_t N>
ACCEL_DECORATION
T Array<T,N>::max() const {
    T res = std::numeric_limits<T>::min();
    for (size_t ii = 0; ii < N; ii++) { res = (res > m_data[ii] ? res : m_data[ii]); }
    return res;
}

template<typename T, size_t N>
ACCEL_DECORATION
T Array<T,N>::min() const {
    T res = std::numeric_limits<T>::max();
    for (size_t ii = 0; ii < N; ii++) { res = (res < m_data[ii] ? res : m_data[ii]); }
    return res;
}

template<typename T, size_t N>
ACCEL_DECORATION
T Array<T,N>::absMax() const {
    T res = 0;
    for (size_t ii = 0; ii < N; ii++) { res = (res > abs(m_data[ii]) ? res : abs(m_data[ii])); }
    return res;
}

template<typename T, size_t N>
ACCEL_DECORATION
T Array<T,N>::sum() const {
    T res = 0;
    for (size_t ii = 0; ii < N; ii++) { res += m_data[ii]; }
    return res;
}

template<typename T, size_t N>
ACCEL_DECORATION
T Array<T,N>::norm() const {
    return std::sqrt(this->dot(*this));
}

template<typename T, size_t N>
ACCEL_DECORATION
T Array<T,N>::product() const {
    T prod = 1.0;
    for (int ii = 0; ii < DIM; ii++) { prod *= m_data[ii]; }
    return prod;
}

template<typename T, size_t N>
ACCEL_DECORATION
T Array<T,N>::dot(const Array<T,N> &a_rhs) const {
    T res = 0;
    for (size_t ii = 0; ii < N; ii++) { res += m_data[ii]*a_rhs[ii]; }
    return res;
}

// Deprecated
// template<typename T, size_t N>
// ACCEL_DECORATION
// void Array<T,N>::reset() { this->fill(0); }

template<typename T, size_t N>
ACCEL_DECORATION
void Array<T,N>::fill(T value) { std::fill_n(m_data,N,value); }

template<typename T, size_t N>
ACCEL_DECORATION
T* Array<T,N>::data() { return m_data; }

template<typename T, size_t N>
ACCEL_DECORATION
const T* Array<T,N>::data() const { return m_data; }

template<typename T, size_t N>
ACCEL_DECORATION
T& Array<T,N>::operator[](size_t i)
{
    PROTO_ASSERT(i < N, "Array::operator[] | Error: index out of bounds.");
    return m_data[i];
}

template<typename T, size_t N>
ACCEL_DECORATION
const T& Array<T,N>::operator[](size_t i) const
{
    PROTO_ASSERT(i < N, "Array::operator[] | Error: index out of bounds.");
    return m_data[i];
}

template<typename T, size_t N>
ACCEL_DECORATION
ArrayIterator<T,N> Array<T,N>::begin() { return ArrayIterator<T,N>(*this, 0); }

template<typename T, size_t N>
ACCEL_DECORATION
ConstArrayIterator<T,N> Array<T,N>::begin() const
{ return ConstArrayIterator<T,N>(*this, 0); }

template<typename T, size_t N>
ACCEL_DECORATION
ArrayIterator<T,N> Array<T,N>::end() { return ArrayIterator<T,N>(*this, N); }

template<typename T, size_t N>
ACCEL_DECORATION
ConstArrayIterator<T,N> Array<T,N>::end() const
{ return ConstArrayIterator<T,N>(*this, N); }

template<typename T, size_t N>
std::string Array<T,N>::str() const {
    std::ostringstream oss;
    oss << "(";
    for (int dir = 0; dir < DIM; dir++)
    {
        oss << m_data[dir];
        if (dir < DIM-1) { oss << ", "; }
    }
    oss << ")";
    return oss.str();
}

template<typename T, size_t N>
void Array<T,N>::print(std::ostream& stream) const {
    stream << str() << std::endl;
}

// Premultiplication by a scalar
template<typename T, size_t N>
ACCEL_DECORATION
Array<T,N>& operator*(int a_scale, Array<T,N>& a_array) {
    for (size_t ii = 0; ii < N; ii++)
    {
        a_array[ii] *= a_scale;
    }
    return a_array;
}

template<typename T, size_t N>
ACCEL_DECORATION
Array<T,N>& operator*(double a_scale, Array<T,N>& a_array) {
    for (size_t ii = 0; ii < N; ii++)
    {
        a_array[ii] *= a_scale;
    }
    return a_array;
}

// Unary negation
template<typename T, size_t N>
ACCEL_DECORATION
Array<T,N>& operator-(Array<T,N>& a_array) {
    for (size_t ii = 0; ii < N; ii++) { a_array[ii] = -a_array[ii]; }
    return a_array;
}

// Stream Operator
template<typename T, size_t N>
std::ostream& operator<<(std::ostream& a_os, const Array<T,N>& a_array)
{
    a_os << a_array.str();
    return a_os;
}
} //end namespace Proto
