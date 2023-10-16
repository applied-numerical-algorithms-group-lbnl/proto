template<typename T, unsigned int C, MemType MEM, Centering CTR>
void PeriodicBC<T,C,MEM,CTR>::apply(LevelBoxData<T,C,MEM,CTR>& a_data) const
{
    a_data.exchange();
}

