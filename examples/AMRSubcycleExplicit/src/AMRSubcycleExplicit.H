template
<class OpType, typename T,unsigned int C,Memtype MEM>
class AMRRKExplicit
{
public:
  inline void define
        (AMRRKExplicit<OpType,T,C,MEM>& a_coarser,
         Point& a_refRatio,
         int timeRefRatio);

  inline void advanceLevel;
  inline void interpCoarse;
  inline void subcycle(double& a_dt);
  
protected:
  OpType m_op;
  shared_ptr<AMRData<T,C,MEM> > m_dataPtr;
  shared_ptr<AMRRKExplicit<OpType,T,C,MEM> > m_finer;
  LevelFluxRegister<T,C,MEM> m_register;
  LevelBoxData<T,RKORDER*C,MEM> m_UTaylor;
  LevelBoxData<T,RKORDER*C,MEM> m_UTaylorCoarsened;
  vector<vector<Box> > m_adjacentBoxes;
  InterpStencil<double> m_cfInterp;
  Point m_ghost;
  int m_timeRefRatio;
  Point m_refRatio;
  int m_numSteps = 0;
  int m_regridInterval = 10000;
  constexpr array<double,RKSTAGES> m_stageWgt = {1./6.,1./3.,1./3.,1./6.};
  constexpr array<double,RKSTAGES> m_stageTimes = {0.,1./2.,1./2.,1.0};
}
