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
    array<double,RKSTAGES>      m_stageWgt = {1./6.,1./3.,1./3.,1./6.};
    array<double,RKSTAGES+1>    m_stageTimes = {0.,1./2.,1./2.,1.,1.};
    
    // data defined on coarse levels
    Point   m_spaceRefRatio;
    int     m_timeRefRatio;
    LevelBoxData<T, RKORDER*C, MEM>             m_UTaylor;
    std::shared_ptr<AMRRK4<OpType, T, C, MEM>>  m_fine;

    // data defined on fine levels
    Point               m_spaceRefRatioCrse;
    int                 m_timeRefRatioCrse;
    InterpStencilTP<T>  m_cfInterp
    LevelBoxData<T, RKORDER*C, MEM>   m_UTaylorCrse;
};
#include "AMRRK4Implem.H"
#endif // end include guard

