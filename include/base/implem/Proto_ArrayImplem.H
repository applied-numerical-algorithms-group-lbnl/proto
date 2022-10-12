template<typename T, size_t N>
ACCEL_DECORATION
Array<T,N>::Array() {
    std::fill_n(m_data,T(0),N);
}
template<typename T, size_t N>
ACCEL_DECORATION
Array<T,N>::Array(const Array<T,N> &arr) {
    proto_memcpy(arr.data(),m_data,N*sizeof(T));
}
template<typename T, size_t N>
ACCEL_DECORATION
Array<T,N>& Array<T,N>::operator=(const Array<T,N> &arr) {
    proto_memcpy(arr.data(),m_data,N*sizeof(T));
}
template<typename T, size_t N>
ACCEL_DECORATION
Array<T,N>::Array(const Array<T,N> &&arr) {
    proto_memcpy(arr.data(),m_data,N*sizeof(T));
    arr.reset();
}
template<typename T, size_t N>
ACCEL_DECORATION
Array<T,N>& Array<T,N>::operator=(const Array<T,N> &&arr) {
    proto_memcpy(arr.data(),m_data,N*sizeof(T));
    arr.reset();
}
template<typename T, size_t N>
Array<T,N>& Array<T,N>::operator=(const std::array<T,N> &arr) {
    proto_memcpy<HOST>(arr.data(),m_data,N*sizeof(T));
}
template<typename T, size_t N>
Array<T,N>::Array(const std::array<T,N> &arr) {
    proto_memcpy<HOST>(arr.data(),m_data,N*sizeof(T));
}
template<typename T, size_t N>
Array<T,N>::Array(const Point &a_pt) {
    for (size_t i=0; i<min(DIM,N); i++)
        (*this)[i] = a_pt[i];
}
template<typename T, size_t N>
ACCEL_DECORATION
inline Array<T,N> Array<T,N>::operator+(const Array<T,N> &a_rhs) const {
    Array<T,N> ret = *this;
    for (size_t i=0; i<N; i++)
        ret[i] += a_rhs[i];
    return ret;
}
template<typename T, size_t N>
ACCEL_DECORATION
inline Array<T,N> Array<T,N>::operator+(T val) const {
    Array<T,N> ret = *this;
    for (auto elt : ret)
        elt += val;
    return ret;
}
template<typename T, size_t N>
ACCEL_DECORATION
inline void Array<T,N>::operator+=(const Array<T,N> &a_rhs) {
    for (size_t i=0; i<N; i++)
        (*this)[i] += a_rhs[i];
}
template<typename T, size_t N>
ACCEL_DECORATION
inline void Array<T,N>::operator+=(T val) {
    for (auto elt : *this)
        elt += val;
}

template<typename T, size_t N>
ACCEL_DECORATION
inline Array<T,N> Array<T,N>::operator-(const Array<T,N> &a_rhs) const {
    Array<T,N> ret = *this;
    for (size_t i=0; i<N; i++)
        ret[i] -= a_rhs[i];
    return ret;
}
template<typename T, size_t N>
ACCEL_DECORATION
inline Array<T,N> Array<T,N>::operator-(T val) const {
    Array<T,N> ret = *this;
    for (auto elt : ret)
        elt -= val;
    return ret;
}
template<typename T, size_t N>
ACCEL_DECORATION
inline void Array<T,N>::operator-=(const Array<T,N> &a_rhs) {
    for (size_t i=0; i<N; i++)
        (*this)[i] -= a_rhs[i];
}
template<typename T, size_t N>
ACCEL_DECORATION
inline void Array<T,N>::operator-=(T val) {
    for (auto elt : *this)
        elt -= val;
}

template<typename T, size_t N>
ACCEL_DECORATION
inline Array<T,N> Array<T,N>::operator*(const Array<T,N> &a_rhs) const {
    Array<T,N> ret = *this;
    for (size_t i=0; i<N; i++)
        ret[i] *= a_rhs[i];
    return ret;
}
template<typename T, size_t N>
ACCEL_DECORATION
inline Array<T,N> Array<T,N>::operator*(T val) const {
    Array<T,N> ret = *this;
    for (auto elt : ret)
        elt *= val;
    return ret;
}
template<typename T, size_t N>
ACCEL_DECORATION
inline void Array<T,N>::operator*=(const Array<T,N> &a_rhs) {
    for (size_t i=0; i<N; i++)
        (*this)[i] *= a_rhs[i];
}
template<typename T, size_t N>
ACCEL_DECORATION
inline void Array<T,N>::operator*=(T val) {
    for (auto elt : *this)
        elt *= val;
}

// check for division by 0
template<typename T, size_t N>
ACCEL_DECORATION
inline Array<T,N> Array<T,N>::operator/(const Array<T,N> &a_rhs) const {
    Array<T,N> ret = *this;
    for (size_t i=0; i<N; i++) {
        PR_assert(a_rhs[i]!=0)
        ret[i] /= a_rhs[i];
    }
    return ret;
}
template<typename T, size_t N>
ACCEL_DECORATION
inline Array<T,N> Array<T,N>::operator/(T val) const {
    Array<T,N> ret = *this;
    for (auto elt : ret) {
        PR_assert(val!=0)
        elt /= val;
    }
    return ret;
}
template<typename T, size_t N>
ACCEL_DECORATION
inline void Array<T,N>::operator/=(const Array<T,N> &a_rhs) {
    for (size_t i=0; i<N; i++) {
        PR_assert(a_rhs[i]!=0)
        (*this)[i] /= a_rhs[i];
    }
}
template<typename T, size_t N>
ACCEL_DECORATION
inline void Array<T,N>::operator/=(T val) {
    for (auto elt : *this) {
        PR_assert(val!=0)
        elt /= val;
    }
}
template<typename T, size_t N>
ACCEL_DECORATION
inline T Array<T,N>::max() const {
    Reduction<T,Max> rxn;
    rxn.reduce(m_data);
    return rxn.localFetch();
}
template<typename T, size_t N>
ACCEL_DECORATION
inline T Array<T,N>::min() const {
    Reduction<T,Min> rxn;
    rxn.reduce(m_data);
    return rxn.localFetch();
}
template<typename T, size_t N>
ACCEL_DECORATION
inline T Array<T,N>::absMax() const {
    Reduction<T,Abs> rxn;
    rxn.reduce(m_data);
    return rxn.localFetch();
}
template<typename T, size_t N>
ACCEL_DECORATION
inline T Array<T,N>::sum() const {
    Reduction<T,Sum> rxn;
    rxn.reduce(m_data);
    return rxn.localFetch();
}
template<typename T, size_t N>
ACCEL_DECORATION
inline T Array<T,N>::norm() const {
    T res = 0;
    for (auto elt : *this)
        res += elt*elt;
    return sqrt(res);
}
template<typename T, size_t N>
ACCEL_DECORATION
inline void Array<T,N>::reset() { m_data = nullptr; }
template<typename T, size_t N>
ACCEL_DECORATION
inline bool Array<T,N>::empty() { m_data == nullptr; }
template<typename T, size_t N>
ACCEL_DECORATION
inline void Array<T,N>::fill(T value) { std::fill_n(m_data,N,value); }
template<typename T, size_t N>
ACCEL_DECORATION
inline const T* Array<T,N>::data() { return &(m_data[0]); }
template<typename T, size_t N>
ACCEL_DECORATION
inline const T* Array<T,N>::data() const { return &(m_data[0]); }
template<typename T, size_t N>
ACCEL_DECORATION
inline T& Array<T,N>::operator[](size_t i) { return std::reference_wrapper<T>(m_data[i]); }
template<typename T, size_t N>
ACCEL_DECORATION
inline T& Array<T,N>::operator[](size_t i) const { return std::reference_wrapper<T>(m_data[i]); }
template<typename T, size_t N>
ACCEL_DECORATION
inline T* Array<T,N>::begin() const { return m_data; }
template<typename T, size_t N>
ACCEL_DECORATION
inline const T* Array<T,N>::cbegin() const { return m_data; }
template<typename T, size_t N>
ACCEL_DECORATION
inline T* Array<T,N>::end() const { return m_data+N; }
template<typename T, size_t N>
ACCEL_DECORATION
inline const T* Array<T,N>::cend() const { return m_data+N; }

// Premultiplication by a scalar
template<typename T, unsigned int N>
ACCEL_DECORATION
inline Array<T,N> operator*(T a_scale, Array<T,N> a_array) {
    for (auto elt : a_array)
        elt *= a_scale;
}

// Unary negation
template<typename T, unsigned int N>
ACCEL_DECORATION
inline Array<T,N> operator-(Array<T,N> a_array) {
    for (auto elt : a_array)
        elt = -elt;
}
