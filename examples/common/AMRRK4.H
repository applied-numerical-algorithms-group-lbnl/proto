
#ifndef _PROTO_AMR_RK4__
#define _PROTO_AMR_RK4__
#define RKORDER 4
#define RKSTAGES 4
#include "ProtoAMR.H"

template<class OpType, typename T,unsigned int C=1,MemType MEM=MEMTYPE_DEFAULT>
class AMRRK4
{
    public:

    inline AMRRK4(){}
    inline AMRRK4(AMRData<T, C, MEM>& a_data, std::array<T, DIM> a_dx);
    inline void define(AMRData<T, C, MEM>& a_data, std::array<T, DIM> a_dx);
    inline void advance(int a_level, T a_dt, T a_chi);

    private:
    
    unsigned int m_regridBaseLevel = 0;
    unsigned int m_bufferSize = 2;
    std::array<double,RKSTAGES>      m_stageWeights = {1./6.,1./3.,1./3.,1./6.};
    std::array<double,RKSTAGES+1>    m_stageTimes = {0.,1./2.,1./2.,1.,1.};
    
    AMRData<T, C, MEM>* m_data;
    
    std::vector<int>    m_timeRefRatios;
    std::vector<Point>  m_spaceRefRatios;
    std::vector<OpType> m_ops;
    std::vector<InterpStencilTP<T>> m_cfInterp;
    std::vector<LevelBoxData<T, C, MEM>> m_UStage;
    std::vector<LevelBoxData<T, C, MEM>> m_RHSTotal;
    std::vector<LevelBoxData<T, C*RKORDER, MEM>> m_UTaylor;
    std::vector<LevelBoxData<T, C*RKORDER, MEM>> m_UTaylorC;
    std::vector<LevelFluxRegister<T, C, MEM>> m_fluxRegisters;
}


#include "AMRRK4Implem.H"
#endif //end include guard
