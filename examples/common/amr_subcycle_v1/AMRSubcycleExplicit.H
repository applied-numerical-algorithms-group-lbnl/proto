#ifndef _PROTO_AMRSUBCYCLEEXPLICIT_
#define _PROTO_AMRSUBCYCLEEXPLICIT_
#define RKORDER 4
#define RKSTAGES 4
#include "ProtoAMR.H"

template<class OpType, typename T,unsigned int C=1,MemType MEM=MEMTYPE_DEFAULT>
class AMRSubcycleExplicit
{
    public:
    
    AMRSubcycleExplicit(){};
    AMRSubcycleExplicit(
            std::shared_ptr<AMRData<T,C,MEM>>   a_dataPtr,
            std::shared_ptr<vector<int >>       a_globalStep,
            const T& a_dx,
            int      a_timeRefRatio,
            int      a_level);

    void define(
            shared_ptr<AMRData<T,C,MEM>>    a_dataPtr,
            shared_ptr<vector<int>>         a_globalStep,
            const T& a_dx,
            int      a_timeRefRatio,
            int      a_level);

    void advance(
            LevelFluxRegister<T,C,MEM>& a_coarseLFR,
            const T& a_dt,
            const T& a_chi,
            bool a_willRegridCoarse);
    
    protected:
    
    // data members for any level.
    OpType  m_op;
    T       m_time;
    T       m_dx;
    Point   m_ghost;
    int     m_level;
    int     m_numSteps = 0;
    shared_ptr<AMRData<T,C,MEM> >   m_dataPtr;
    shared_ptr<vector<int> >        m_globalStep;

    array<double,RKSTAGES>      m_stageWgt = {1./6.,1./3.,1./3.,1./6.};
    array<double,RKSTAGES+1>    m_stageTimes = {0.,1./2.,1./2.,1.,1.};

    // Data members that are defined for all but the finest level.
    shared_ptr<AMRSubcycleExplicit<OpType,T,C,MEM> > m_finePtr;
    shared_ptr<LevelFluxRegister<T,C,MEM> > m_register;
    LevelBoxData<T,RKORDER*C,MEM> m_UTaylor; 
    Point m_refRatioFine;
    int m_timeRefRatioFine;
    int m_regridInterval = 10000;

    // Data members that are defined for all but the coarsest level.
    LevelBoxData<T,RKORDER*C,MEM> m_UTaylorCoarsened;
    Point m_refRatioCoarse;
    int m_timeRefRatioCoarse;
    InterpStencil<T> m_cfInterp;
    InterpStencilTP<T> m_cfInterpTP;
};
#include "AMRSubcycleExplicitImplem.H"
#endif
