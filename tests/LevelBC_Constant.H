#pragma once

#include "Proto.H"

template<typename T, unsigned int C, MemType MEM, Centering CTR>
class ConstantBC : public LevelBC<T,C,MEM,CTR>
{
    public:
    
    using LevelBC<T,C,MEM,CTR>::LevelBC;

#if 0
    inline void init()
    {
        m_J.define(this->layout());
        // compute J;
    }
#endif
    inline void apply(LevelBoxData<T,C,MEM,CTR>& a_data) const
    {
#if 0
        this->layout(); // DisjointBoxLayout / ProblemDomain
        this->getLevelData(0); //access external level data
#endif
        a_data.setDomainBoundary(this->getConstant(0));
    }

    // LevelBoxData<double, 1> m_J;
    // LevelBoxData<double, DIM> m_X;
};
