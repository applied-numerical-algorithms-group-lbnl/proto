#ifndef _PROTO_AMR_RK4_
#define _PROTO_AMR_RK4_

#define RKORDER 4
#define RKSTAGES 4
#include "ProtoAMR.H"

using namespace Proto;

template<class OpType, typename T,unsigned int C=1,MemType MEM=MEMTYPE_DEFAULT>
class AMRRK4
{
    public:

    inline AMRRK4(){};
    
    inline AMRRK4(
            AMRData<T,C,MEM>&   a_data,
            std::array<double>  a_dx);

    inline void define(
            AMRData<T,C,MEM>&   a_data,
            std::array<double>  a_dx);

    inline void advance(const T& a_dt);

    private:

    inline AMRRK4(
            AMRData<T,C,MEM>&   a_data,
            std::array<double>  a_dx,
            unsigned int        a_level);
    
    inline void define(
            AMRData<T,C,MEM>&   a_data,
            std::array<double>  a_dx,
            unsigned int        a_level);

    // data defined on all levels
    OpType              m_op;
    AMRData<T, C, MEM>* m_data;
    unsigned int        m_level;
    T                   m_time;
    std::array<T, DIM>  m_dx;

    // data defined on coarse levels
    Point   m_spaceRefRatio;
    int     m_timeRefRatio;
    std::shared_ptr<AMRRK4<OpType, T, C, MEM>>              m_fine;
    std::vector<std::shared_ptr<LevelBoxData<T, C, MEM>>>   m_UTaylor;

    // data defined on fine levels
    Point               m_spaceRefRatioCrse;
    int                 m_timeRefRatioCrse;
    InterpStencilTP<T>  m_cfInterp
    std::vector<std::shared_ptr<LevelBoxData<T, C, MEM>>>   m_UTaylorCrse;
};
#include "AMRRK4Implem.H"
#endif // end include guard

