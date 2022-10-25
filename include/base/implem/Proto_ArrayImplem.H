template<typename T, size_t N>
ACCEL_DECORATION
Array<T,N>::Array() {}

template<typename T, size_t N>
ACCEL_DECORATION
Array<T,N>::Array(T val) {
    for (auto &elt : *this)
        elt = val;
}

template<typename T, size_t N>
ACCEL_DECORATION
Array<T,N>::Array(std::initializer_list<T> list) {
    size_t i=0;
    for (auto &elt : list)
        m_data[i++] = elt;
} 

template<typename T, size_t N>
ACCEL_DECORATION
Array<T,N>::Array(const Array<T,N> &arr) {
    size_t i=0;
    for (auto &elt : *this)
        elt = arr[i++];
}

template<typename T, size_t N>
ACCEL_DECORATION
Array<T,N>& Array<T,N>::operator=(const Array<T,N> &arr) {
    size_t i=0;
    for (auto &elt : *this)
        elt = arr[i++];
    return *this;
}

template<typename T, size_t N>
Array<T,N>& Array<T,N>::operator=(const std::array<T,N> &arr) {
    size_t i=0;
    for (auto &elt : *this)
        elt = arr[i++];
    return *this;
}

template<typename T, size_t N>
Array<T,N>::Array(const std::array<T,N> &arr) {
    size_t i=0;
    for (auto &elt : *this)
        elt = arr[i++];
}

template<typename T, size_t N>
Array<T,N>::Array(const Point &a_pt) {
    for (size_t i=0; i<min(DIM,N); i++)
        m_data[i] = a_pt[i];
}

template<typename T, size_t N>
ACCEL_DECORATION
inline bool Array<T,N>::operator==(const Array<T,N> &a_rhs) const {
    for (size_t i=0; i<N; i++)
        if (m_data[i] != a_rhs[i])
            return false;
    return true;
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
        m_data[i] += a_rhs[i];
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
        m_data[i] -= a_rhs[i];
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
        m_data[i] *= a_rhs[i];
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
        m_data[i] /= a_rhs[i];
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
    T res = std::numeric_limits<T>::min();
    for (auto elt : *this)
        res = (res > elt ? res : elt);
    return res;
}

template<typename T, size_t N>
ACCEL_DECORATION
inline T Array<T,N>::min() const {
    T res = std::numeric_limits<T>::max();
    for (auto elt : *this)
        res = (res < elt ? res : elt);
    return res;
}

template<typename T, size_t N>
ACCEL_DECORATION
inline T Array<T,N>::absMax() const {
    T res = 0;
    for (auto elt : *this)
        res += (elt > 0 ? elt : -elt);
    return res;
}

template<typename T, size_t N>
ACCEL_DECORATION
inline T Array<T,N>::sum() const {
    T res = 0;
    for (auto elt : *this)
        res += elt;
    return res;
}

template<typename T, size_t N>
ACCEL_DECORATION
inline T Array<T,N>::norm() const {
    return sqrt(dot(*this));
}

template<typename T, size_t N>
ACCEL_DECORATION
inline T Array<T,N>::dot(const Array<T,N> &a_rhs) const {
    T res = 0; size_t i = 0;
    for (auto &elt : a_rhs)
        res += m_data[i++]*elt;
    return res;
}

template<typename T, size_t N>
ACCEL_DECORATION
inline void Array<T,N>::reset() { this->fill(0); }

template<typename T, size_t N>
ACCEL_DECORATION
inline void Array<T,N>::fill(T value) { std::fill_n(m_data,N,value); }

template<typename T, size_t N>
ACCEL_DECORATION
inline T* Array<T,N>::data() { return const_cast<T*>(&(m_data[0])); }

template<typename T, size_t N>
ACCEL_DECORATION
inline T* Array<T,N>::data() const { return const_cast<T*>(&(m_data[0])); }

template<typename T, size_t N>
ACCEL_DECORATION
inline T& Array<T,N>::operator[](size_t i) { return const_cast<T&>(m_data[i]); }

template<typename T, size_t N>
ACCEL_DECORATION
inline T& Array<T,N>::operator[](size_t i) const { return const_cast<T&>(m_data[i]); }

template<typename T, size_t N>
ACCEL_DECORATION
inline T* Array<T,N>::begin() const { return this->data(); }

template<typename T, size_t N>
ACCEL_DECORATION
inline const T* Array<T,N>::cbegin() const { return this->data(); }

template<typename T, size_t N>
ACCEL_DECORATION
inline T* Array<T,N>::end() const { return this->data()+N; }

template<typename T, size_t N>
ACCEL_DECORATION
inline const T* Array<T,N>::cend() const { return this->data()+N; }

template<typename T, size_t N>
inline void Array<T,N>::print() const {
    std::cout << "( ";
    for (auto &elt : *this)
        std::cout << elt << ' ';
    std::cout << ')' << std::endl;
}

// Premultiplication by a scalar
template<typename T, size_t N>
ACCEL_DECORATION
inline Array<T,N>& operator*(T a_scale, Array<T,N> &a_array) {
    for (auto elt : a_array)
        elt *= a_scale;
    return a_array;
}

// Unary negation
template<typename T, size_t N>
ACCEL_DECORATION
inline Array<T,N>& operator-(Array<T,N> &a_array) {
    for (auto elt : a_array)
        elt = -elt;
    return a_array;
}
