
#ifndef _PROTO_AMR_RK4__
#define _PROTO_AMR_RK4__
#define RKORDER 4
#define RKSTAGES 4
#include "Proto.H"

extern int TIME_STEP;
extern int PLOT_NUM;

template<class OpType, typename T,unsigned int C=1,MemType MEM=MEMTYPE_DEFAULT>
class AMRRK4
{
    public:

    inline AMRRK4(){}
    inline AMRRK4(AMRData<T, C, MEM>& a_data, const Array<T, DIM>& a_dx);
    inline void define(AMRData<T, C, MEM>& a_data, const Array<T, DIM>& a_dx);
    inline void advance(T a_dt);

    private:
   
    inline void define(AMRData<T, C, MEM>& a_data, const Array<T, DIM>& a_dx, int a_level);
    inline void advance(int a_level, T a_dt, T a_chi);
    inline void interpBoundaries(int a_level, int a_stage, T a_chi);
    inline void stageUpdate(int a_level, int a_stage, T a_dt);
    inline void regrid(int a_baseLevel);

    unsigned int m_regridBaseLevel = 0;
    unsigned int m_bufferSize = 2;
    Array<double,RKSTAGES>      m_stageWeights = {1./6.,1./3.,1./3.,1./6.};
    Array<double,RKSTAGES+1>    m_stageTimes = {0.,1./2.,1./2.,1.,1.};
    
    AMRData<T, C, MEM>* m_data;
    
    std::vector<int>    m_timeRefRatios;
    std::vector<Point>  m_spaceRefRatios;
    std::vector<OpType> m_ops;
    std::vector<InterpStencilTP<T>> m_cfInterp;
    std::vector<std::shared_ptr<LevelBoxData<T, C, MEM>>> m_UStage;
    std::vector<std::shared_ptr<LevelBoxData<T, C, MEM>>> m_RHSTotal;
    std::vector<std::shared_ptr<LevelBoxData<T, C*RKORDER, MEM>>> m_UTaylor;
    std::vector<std::shared_ptr<LevelBoxData<T, C*RKORDER, MEM>>> m_UTaylorC;
    std::vector<std::shared_ptr<LevelFluxRegister<T, C, MEM>>> m_fluxRegisters;
};


#include "AMRRK4Implem.H"
#endif //end include guard
