template<typename T, size_t N, MemType Mem>
Array<T,N,Mem>::Array() {
    m_data = (T*)proto_malloc<Mem>(N*sizeof(T));
}
template<typename T, size_t N, MemType Mem>
template<MemType Ory>
Array<T,N,Mem>::Array(const Array<T,N,Ory> &arr) : Array<T,N,Mem>() {
    proto_memcpy<Ory,Mem>(arr.data(),m_data,N*sizeof(T));
}
template<typename T, size_t N, MemType Mem>
template<MemType Ory>
Array<T,N,Mem>& Array<T,N,Mem>::operator=(const Array<T,N,Ory> &arr) {
    if (m_data == nullptr) {
        proto_malloc<Mem>(m_data,N*sizeof(T));
    }
    proto_memcpy<Ory,Mem>(arr.data(),m_data,N*sizeof(T));
}
template<typename T, size_t N, MemType Mem>
template<MemType Ory>
Array<T,N,Mem>::Array(const Array<T,N,Ory> &&arr) : Array<T,N,Mem>() {
    proto_memcpy<Ory,Mem>(arr.data(),m_data,N*sizeof(T));
    arr.reset();
}
template<typename T, size_t N, MemType Mem>
template<MemType Ory>
Array<T,N,Mem>& Array<T,N,Mem>::operator=(const Array<T,N,Ory> &&arr) {
    if (m_data == nullptr) {
        proto_malloc<Mem>(m_data,N*sizeof(T));
    }
    proto_memcpy<Ory,Mem>(arr.data(),m_data,N*sizeof(T));
    arr.reset();
}
template<typename T, size_t N, MemType Mem>
Array<T,N,Mem>& Array<T,N,Mem>::operator=(const std::array<T,N> &arr) {
    if (m_data == nullptr) {
        proto_malloc<Mem>(m_data,N*sizeof(T));
    }
    proto_memcpy<HOST,Mem>(arr.data(),m_data,N*sizeof(T));
}
template<typename T, size_t N, MemType Mem>
Array<T,N,Mem>::Array(const std::array<T,N> &arr) : Array<T,N,Mem>() {
    proto_memcpy<HOST,Mem>(arr.data(),m_data,N*sizeof(T));
}
template<typename T, size_t N, MemType Mem>
Array<T,N,Mem>::~Array() { 
    proto_free<Mem>(m_data); 
}
template<typename T, size_t N, MemType Mem>
ACCEL_DECORATION
inline Array<T,N,Mem> operator+(const Array<T,N,Mem>& a_rhs) const {
    Array<T,N,Mem> ret;
    for (size_t i=0; i<N; i++)
        ret[i] = (*this)[i] + a_rhs[i];
    return ret;
}
template<typename T, size_t N, MemType Mem>
ACCEL_DECORATION
inline void Array<T,N,Mem>::reset() { m_data = nullptr; }
template<typename T, size_t N, MemType Mem>
ACCEL_DECORATION
inline void Array<T,N,Mem>::fill(T value) { std::fill_n(this->data(),N,value); }
template<typename T, size_t N, MemType Mem>
ACCEL_DECORATION
inline T* Array<T,N,Mem>::data() { return m_data; }
template<typename T, size_t N, MemType Mem>
ACCEL_DECORATION
inline T* Array<T,N,Mem>::data() const { return m_data; }
template<typename T, size_t N, MemType Mem>
ACCEL_DECORATION
inline T& Array<T,N,Mem>::operator[](size_t i) { return std::reference_wrapper<T>(m_data[i]); }
template<typename T, size_t N, MemType Mem>
ACCEL_DECORATION
inline T& Array<T,N,Mem>::operator[](size_t i) const { return std::reference_wrapper<T>(m_data[i]); }
