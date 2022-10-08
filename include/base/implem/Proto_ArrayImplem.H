template<typename T, size_t N, MemType Mem>
ACCEL_DECORATION
Array<T,N,Mem>::Array() {
    m_data = new T[N];
}
template<typename T, size_t N, MemType Mem>
ACCEL_DECORATION
Array<T,N,Mem>::Array(const Array<T,N,Mem> &arr) : Array<T,N,Mem>() {
    memcpy(m_data,arr.data(),N*sizeof(T));
}
template<typename T, size_t N, MemType Mem>
ACCEL_DECORATION
Array<T,N,Mem>& Array<T,N,Mem>::operator=(const Array<T,N,Mem> &arr) {
    memcpy(m_data,arr.data(),N*sizeof(T));
}
template<typename T, size_t N, MemType Mem>
Array<T,N,Mem>& Array<T,N,Mem>::operator=(const std::array<T,N> &arr) {
    if (Mem==DEVICE) {
        protoMemcpy(Mem,m_data,arr.data(),N*sizeof(T),protoMemcpyHostToDevice);
    } else {
        memcpy(m_data,arr.data(),N*sizeof(T));
    }
}
template<typename T, size_t N, MemType Mem>
Array<T,N,Mem>::Array(const std::array<T,N> &arr) : Array<T,N,Mem>() {
    if (Mem==DEVICE) {
        protoMemcpy(Mem,m_data,arr.data(),N*sizeof(T),protoMemcpyHostToDevice);
    } else {
        memcpy(m_data,arr.data(),N*sizeof(T));
    }
}
template<typename T, size_t N, MemType Mem>
ACCEL_DECORATION
Array<T,N,Mem>::~Array() { delete m_data; }
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
