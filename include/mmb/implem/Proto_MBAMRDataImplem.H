template<typename T, unsigned int C, MemType MEM, Centering CTR>
MBAMRData<T,C,MEM,CTR>::MBAMRData(){}
template<typename T, unsigned int C, MemType MEM, Centering CTR>
MBAMRData<T,C,MEM,CTR>::MBAMRData(
        const MBAMRGrid& a_grid,
        Array<Point, DIM+1> a_ghost,
        Point a_boundGhost = Point::Zeros()){}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void MBAMRData<T,C,MEM,CTR>::define(
        const MBAMRGrid& a_grid,
        Array<Point, DIM+1> a_ghost,
        Point a_boundGhost = Point::Zeros()){}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void MBAMRData<T,C,MEM,CTR>::setVal(T a_value){}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void MBAMRData<T,C,MEM,CTR>::setRandom(T a_low, T a_high){}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void MBAMRData<T,C,MEM,CTR>::exchange(){}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
MBLevelBoxData<T,C,MEM,CTR>&
MBAMRData<T,C,MEM,CTR>::getLevel(unsigned int a_level){}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
const MBLevelBoxData<T,C,MEM,CTR>&
MBAMRData<T,C,MEM,CTR>::getLevel(unsigned int a_level) const{}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
AMRData<T,C,MEM,CTR>&
MBAMRData<T,C,MEM,CTR>::getBlock(unsigned int a_block){}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
const AMRData<T,C,MEM,CTR>&
MBAMRData<T,C,MEM,CTR>::getBlock(unsigned int a_block) const{}

// Operators
template<typename T, unsigned int C, MemType MEM, Centering CTR>
void MBAMRData<T,C,MEM,CTR>::averageDown(){}        

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void MBAMRData<T,C,MEM,CTR>::increment(
        const MBAMRData<T,C,MEM,CTR>& a_data, T a_scale){}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void MBAMRData<T,C,MEM,CTR>::operator*=(T a_value){}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void MBAMRData<T,C,MEM,CTR>::operator+=(T a_value){}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
T MBAMRData<T,C,MEM,CTR>::integrate(Array<T,DIM> a_cdx, unsigned int a_c = 0){}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
T MBAMRData<T,C,MEM,CTR>::integrate(T a_cdx, unsigned int a_c = 0){}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
T MBAMRData<T,C,MEM,CTR>::integrateAbs(Array<T,DIM> a_cdx, unsigned int a_c = 0){}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
T MBAMRData<T,C,MEM,CTR>::integrateAbs(T a_cdx, unsigned int a_c = 0){}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
T MBAMRData<T,C,MEM,CTR>::absMax(unsigned int a_c = 0){}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void MBAMRData<T,C,MEM,CTR>::copyTo(MBAMRData& a_rhs){}
