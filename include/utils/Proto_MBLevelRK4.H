#pragma once
#ifndef _PROTO_MB_LEVEL_RK4_H_
#define _PROTO_MB_LEVEL_RK4_H_

#include "Proto.H"
#include "Proto_CubedSphereShell.H"
//#include "Inputs_Parsing.H"

namespace Proto {

///  Generic Explicit EulerStep Algorithm
/**  Given y'=f(y,t), it computes y_{t+1} from y_t using the following equations:
     k1=dt*f(y_t,t)
     y_{t+1}=y_t+k1
*/

///  Generic Explicit RK4 Algorithm
/**  Given y'=f(y,t), it computes y_{t+1} from y_t using the following equations:
     k1=dt*f(y_t,t)
     k2=dt*f(y_t+0.5*k1,t+0.5*dt)
     k3=dt*f(y_t+0.5*k2,t+0.5*dt)
     k4=dt*f(y_t+k3,t+dt)
     y_{t+1}=y_t+k1/6+k2/3+k3/3+k4/6

*/
///  Generic Explicit RK3 Algorithm
/**  Given y'=f(y,t), it computes y_{t+1} from y_t using the following equations:
     k1=dt*f(y_t,t)
     k2=dt*f(y_t+0.5*k1,t+0.5*dt)
     k3=dt*f(y_t-k1+2*k2,t+dt)
     y_{t+1}=y_t+k1/6+2*k2/3+k3/6

*/

template <
    template<typename, typename, MemType> class OPType,
    class MAP,
    typename T,
    template<typename, unsigned int, typename, MemType, Centering> class BCType = NullBC,
    MemType MEM = MEMTYPE_DEFAULT>
class MBLevelRK4
{
public:
    typedef OPType<T, MAP, MEM> OP; 
    typedef BCType<T,OP::numState(), MAP, MEM, PR_CELL> BC;
    typedef MBLevelOp<OPType, MAP, T , BCType, MEM> LOP;
    typedef MBLevelBoxData<T, OP::numState(), MEM, PR_CELL> LevelStateData;
    typedef MBLevelBoxData<T, 1, MEM, PR_CELL> LevelScalarData;
    
    MBLevelRK4(
            const MBLevelMap<MAP>&      a_map,
            MBInterpOp&                 a_interp);

    MBLevelRK4(const MBLevelMap<MAP>& a_map, const MBInterpOp& a_interp);

    //Compute y_{t+1} at time a_time and with y_t=a_state, and place it in a_state.
  inline void advance(LevelStateData& a_state,
                      LevelScalarData& a_dVolrLev, 
                      T a_dt,
                      T a_time,
                      int a_order);
  
    inline void advance(LevelStateData& a_state, double& a_dt, double a_time = 0.);
    inline LOP& op() { return m_f; }
    inline BC& bc() {return *m_bc; }
     template<unsigned int C>
     inline Array<T,C> getCons()
  {
    Array<T,C> consOut;
    for (int comp = 0; comp < C ;comp++)
      {
        consOut[comp] = m_consSum[comp].fetch();
      }
    return consOut;
  }
protected:
    LevelStateData      m_k;
    LevelStateData      m_k1;
    LevelStateData      m_k2;
    LevelStateData      m_delta;
    LOP                 m_f;
    std::shared_ptr<BC> m_bc;
    LevelStateData      m_U0;
    bool                m_U0Defined = false;
    LevelStateData      m_stageState;
    bool m_stageDefined = false;
    const MBInterpOp&   m_interp;
    MBInterpOp*         m_stageInterp = nullptr;
  MBBoundaryRegister<T,OP::numState(),MEM> m_blockreg;
  Array<Reduction<T,Operation::Sum,MEM>,OP::numState()> m_consSum;
  int m_counter = 0;
};
template <
    template<typename, typename, MemType> class OPType,
    class MAP,
    typename T,
    template<typename, unsigned int, typename, MemType, Centering> class BCType,
    MemType MEM>
MBLevelRK4<OPType, MAP, T, BCType, MEM>::MBLevelRK4(
    const MBLevelMap<MAP>&      a_map,
    const MBInterpOp&  a_interp) : m_interp(a_interp)
{
    Array<Point,DIM+1> ghost;
    ghost.fill(Point::Zeros());
    m_f.define(a_map);
    MBDisjointBoxLayout layout = a_map.layout();
    m_blockreg.define(layout,1,Point::Zeros());
    m_delta.define(a_map.layout(), ghost);
    m_k.define(a_map.layout(), ghost);
    m_k1.define(a_map.layout(), ghost);
    m_k2.define(a_map.layout(), ghost);
    // The MHD advance overload does not use the generic RK4 scratch state.
    m_bc = std::make_shared<BC>(a_map);
    for (int comp = 0; comp < OP::numState();comp++) m_consSum[comp].init();
    m_counter = 0;
}

template <
    template<typename, typename, MemType> class OPType,
    class MAP, typename T,
    template<typename, unsigned int, typename, MemType, Centering> class BCType,
    MemType MEM>
MBLevelRK4<OPType, MAP, T, BCType, MEM>::MBLevelRK4(
    const MBLevelMap<MAP>& a_map, MBInterpOp& a_interp)
    : MBLevelRK4(a_map, static_cast<const MBInterpOp&>(a_interp))
{
    m_stageInterp = &a_interp;
}
template <
    template<typename, typename, MemType> class OPType,
    class MAP,
    typename T,
    template<typename, unsigned int, typename, MemType, Centering> class BCType,
    MemType MEM>
void MBLevelRK4<OPType, MAP, T, BCType, MEM>::
advance(
        LevelStateData& a_state,
        LevelScalarData& a_dVolrLev, 
        T a_dt,
        T a_time,
        int a_order)
{
  PR_TIMERS("MBLevelRK4::advance");
    PROTO_ASSERT(m_stageInterp, "The MHD integrator requires a mutable interpolation operator");
    PROTO_ASSERT(a_order == 1 || a_order == 3 || a_order == 4,
        "Supported MHD temporal orders are 1, 3 and 4");
    if (!m_stageDefined)
    {
        Array<Point, DIM+1> ghost;
        ghost.fill(OP::ghost() + Point::Ones(2));
        ghost[0] = OP::ghost();
        m_stageState.define(a_state.layout(), ghost);
        m_stageDefined = true;
    }
    T sixth = 1, third=1, twoThirds = 2, half = 1;
    sixth/=6; third/=3; twoThirds/=3; half/=2;
    int kstage;
    m_delta.setVal(0.);
    m_blockreg.clear();
    
    bool flagi = true; 
    
    if (a_order == 1)
    {
      T weight = 1.0;
      a_state.copyTo(m_stageState);
      kstage = 0;
        {
          OP::preStageLevel(m_stageState, m_k, a_dVolrLev, *m_stageInterp, m_blockreg, m_consSum,
                        a_time, a_dt, weight, kstage, true);
          m_delta.increment(m_k, weight*a_dt);
          if (flagi)
            {
              OP::StageLevel(m_k,m_stageState, m_blockreg, m_consSum,
                             a_time, a_dt, weight, kstage);
              m_delta.increment(m_k, weight*a_dt);
            }
        }
    }

    if (a_order == 3)
    {
        a_state.copyTo(m_stageState);
        kstage = 0;
        T weight = sixth;

        // Taking the Euler step for artificial viscosity
        OP::preStageLevel(m_stageState, m_k, a_dVolrLev, *m_stageInterp, m_blockreg, m_consSum,
                          a_time, a_dt, weight, kstage, true);
        a_state.increment(m_k, a_dt);

        // Now doing the RK3 steps
        a_state.copyTo(m_stageState);
        OP::preStageLevel(m_stageState, m_k, a_dVolrLev, *m_stageInterp, m_blockreg, m_consSum,
                          a_time, a_dt, weight, kstage, false);
        if (flagi)
          {
            OP::StageLevel
              (m_k, m_stageState, m_blockreg, m_consSum,
               a_time, a_dt, weight, kstage);
            {
              PR_TIMERS("block_increment_0");
              m_delta.increment(m_k, a_dt*weight);
              m_k.copyTo(m_k1);
            }
          }

        a_state.copyTo(m_stageState);
        m_stageState.increment(m_k1, half*a_dt);
        kstage = 1;
        weight = twoThirds;
        OP::preStageLevel
          (m_stageState, m_k, a_dVolrLev, *m_stageInterp, m_blockreg, m_consSum,
           a_time+0.5*a_dt, a_dt, weight, kstage, false);
        if (flagi)
          {
            OP::StageLevel
              (m_k, m_stageState, m_blockreg, m_consSum,
               a_time, a_dt, weight, kstage);
            {
              PR_TIMERS("block_increment_1");
              m_delta.increment(m_k, a_dt*weight);
              m_k.copyTo(m_k2);
            }
          }

        a_state.copyTo(m_stageState);
        m_stageState.increment(m_k1, -a_dt);
        m_stageState.increment(m_k2, 2.0*a_dt);
        kstage = 2;
        weight = sixth;
        OP::preStageLevel(m_stageState, m_k, a_dVolrLev, *m_stageInterp,m_blockreg, m_consSum,
                          a_time+a_dt,a_dt, weight, kstage, false);
        if (flagi)
          {
            OP::StageLevel
              (m_k, m_stageState, m_blockreg, m_consSum,
               a_time, a_dt, weight, kstage);
           {
              PR_TIMERS("block_increment_2");
              m_delta.increment(m_k, a_dt*weight);
            }
          }
    }

    if (a_order == 4)    
    {
        a_state.copyTo(m_stageState);
        kstage = 0;
        T weight = sixth;
        
        // Taking the Euler step for artificial viscosity
        OP::preStageLevel(m_stageState, m_k, a_dVolrLev, *m_stageInterp, m_blockreg, m_consSum,
                          a_time, a_dt, weight, kstage, true);
        a_state.increment(m_k, a_dt);

        // Now doing the RK4 steps
        a_state.copyTo(m_stageState);
        OP::preStageLevel(m_stageState, m_k, a_dVolrLev, *m_stageInterp, m_blockreg, m_consSum,
                          a_time, a_dt, weight, kstage, false);
        if (flagi)
          {
            OP::StageLevel
              (m_k, m_stageState, m_blockreg, m_consSum,
               a_time, a_dt, weight, kstage);
            {
              PR_TIMERS("block_increment_0");
              m_delta.increment(m_k, a_dt*weight);
            }
          }

        a_state.copyTo(m_stageState);
        m_stageState.increment(m_k, half*a_dt);
        kstage = 1;
        weight = third;
        OP::preStageLevel
          (m_stageState, m_k, a_dVolrLev, *m_stageInterp, m_blockreg, m_consSum,
           a_time+0.5*a_dt, a_dt, weight, kstage, false);
        if (flagi)
          {
            OP::StageLevel
              (m_k, m_stageState, m_blockreg, m_consSum,
               a_time, a_dt, weight, kstage);
            {
              PR_TIMERS("block_increment_1");
              m_delta.increment(m_k, a_dt*weight);
            }
          }

        a_state.copyTo(m_stageState);
        m_stageState.increment(m_k, half*a_dt);
        kstage = 2;
        weight = third;
        OP::preStageLevel(m_stageState, m_k, a_dVolrLev, *m_stageInterp,m_blockreg, m_consSum,
                          a_time+0.5*a_dt,a_dt, weight, kstage, false);
        if (flagi)
          {
            OP::StageLevel
              (m_k, m_stageState, m_blockreg, m_consSum,
               a_time, a_dt, weight, kstage);
           {
              PR_TIMERS("block_increment_2");
              m_delta.increment(m_k, a_dt*weight);
            }
          }

        a_state.copyTo(m_stageState);
        m_stageState.increment(m_k,a_dt);
        kstage = 3;
        weight = sixth;
        OP::preStageLevel(m_stageState, m_k, a_dVolrLev, *m_stageInterp, m_blockreg, m_consSum,
                          a_time+a_dt,a_dt, weight, kstage, false);
        if (flagi)
          {
            OP::StageLevel
              (m_k, m_stageState, m_blockreg, m_consSum, a_time, a_dt, weight, kstage);
            {
              PR_TIMERS("block_increment_3");
              m_delta.increment(m_k, a_dt*weight);
            }
          }
    }

    // Update State.
    
    //int boundary_reflux = ParseInputs::get_boundary_reflux();
    int boundary_reflux = 1;
    if (boundary_reflux == 1){
      {
        PR_TIMERS("block_reflux");
        m_blockreg.reflux(m_delta);
      }
    }
    a_state.increment(m_delta, 1.0);
    m_counter++;
}
template <
    template<typename, typename, MemType> class OPType,
    class MAP,
    typename T,
    template<typename, unsigned int, typename, MemType, Centering> class BCType,
    MemType MEM>
void MBLevelRK4<OPType, MAP, T, BCType, MEM>::advance(
    LevelStateData& a_state, 
    double&         a_dt,
    double          a_time)
{
    if (!m_U0Defined)
    {
        Array<Point,DIM+1> ghost;
        ghost.fill(OP::ghost());
        m_U0.define(a_state.layout(), ghost);
        m_U0Defined = true;
    }
    double sixth = 1, third=1, half = 1;
    sixth/=6; third/=3; half/=2;

    m_delta.setVal(0);

    m_f.setTime(a_time);
    m_bc->setTime(a_time);
    m_bc->setRKStage(0);
    m_bc->apply(a_state);
    a_state.exchange();
    m_interp.apply(a_state, a_state);
    m_f(m_k, a_state, a_dt);    // compute k1*dt
    m_delta.increment(m_k, sixth);
    
    m_f.setTime(a_time+half*a_dt);
    m_bc->setTime(a_time + half*a_dt);
    m_bc->setRKStage(1);
    a_state.copyTo(m_U0);
    m_U0.increment(m_k, half);
    m_bc->apply(m_U0);
    m_U0.exchange();
    m_interp.apply(m_U0, m_U0);
    m_f(m_k, m_U0, a_dt);  // compute k2*dt
    m_delta.increment(m_k, third);
    
    m_bc->setRKStage(2);
    a_state.copyTo(m_U0);
    m_U0.increment(m_k, half);
    m_bc->apply(m_U0);
    m_U0.exchange();
    m_interp.apply(m_U0, m_U0);
    m_f(m_k, m_U0, a_dt);  // compute k3*dt
    m_delta.increment(m_k, third);
    
    m_f.setTime(a_time+a_dt);
    m_bc->setTime(a_time + a_dt);
    m_bc->setRKStage(3);
    a_state.copyTo(m_U0);
    m_U0.increment(m_k);
    m_bc->apply(m_U0);
    m_U0.exchange();
    m_interp.apply(m_U0, m_U0);
    m_f(m_k, m_U0, a_dt); // compute k4*dt
    m_delta.increment(m_k, sixth);
    
    a_state.increment(m_delta);
}

} //end Proto namespace

#endif //end include guard
