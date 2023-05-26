#pragma once

#include "Proto.H"

template<typename T, unsigned int C, MemType MEM, Centering CTR>
class ConstantBC : public LevelBC<T,C,MEM,CTR>
{
    public:
    
    using LevelBC<T,C,MEM,CTR>::LevelBC;

    inline void apply(LevelBoxData<T,C,MEM,CTR>& a_data) const
    {
        a_data.setBoundary(this->getConstant(0));
    }
};
