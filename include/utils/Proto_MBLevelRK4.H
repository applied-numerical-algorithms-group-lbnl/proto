#pragma once
#ifndef _PROTO_MB_LEVEL_RK4_H_
#define _PROTO_MB_LEVEL_RK4_H_

#include "Proto.H"

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

template <
    template<typename, template<MemType> typename, MemType> class OPType,
    template<MemType> class MAP,
    typename T,
    template<typename, unsigned int, template<MemType> typename, MemType, Centering> class BCType = NullBC,
    MemType MEM = MEMTYPE_DEFAULT>
class MBLevelRK4
{
public:
    typedef OPType<T, MAP, MEM> OP; 
    typedef BCType<T,OP::numState(), MAP, MEM, PR_CELL> BC;
    typedef MBLevelOp<OPType, MAP, double, BCType, MEM> LOP;
    typedef MBLevelBoxData<T, OP::numState(), MEM, PR_CELL> LevelStateData;
    typedef MBLevelBoxData<T, 1, MEM, PR_CELL> LevelScalarData;
    
    MBLevelRK4(
            const MBLevelMap<MAP>&      a_map,
            MBInterpOp&           a_interp);

    //Compute y_{t+1} at time a_time and with y_t=a_state, and place it in a_state.
    inline void advance(LevelStateData& a_state,LevelScalarData& a_dVolrLev, double& a_dt, double a_time = 0., int a_order = 4);
    inline LOP& op() { return m_f; }
    inline BC& bc() {return *m_bc; }
protected:
    LevelStateData      m_k;
    LevelStateData      m_delta;
    LOP                 m_f;
    std::shared_ptr<BC> m_bc;
    LevelStateData      m_U0;
    MBInterpOp&   m_interp;
};
template <
    template<typename, template<MemType> typename, MemType> class OPType,
    template<MemType> class MAP,
    typename T,
    template<typename, unsigned int, template<MemType> typename, MemType, Centering> class BCType,
    MemType MEM>
MBLevelRK4<OPType, MAP, T, BCType, MEM>::MBLevelRK4(
    const MBLevelMap<MAP>&      a_map,
    MBInterpOp&           a_interp) : m_interp(a_interp)
{
    Array<Point,DIM+1> ghost;
    ghost.fill(Point::Zeros());
    m_f.define(a_map);
    m_delta.define(a_map.layout(), ghost);
    m_k.define(a_map.layout(), ghost);
    ghost.fill(OP::ghost() + Point::Ones(2));
    ghost[0] = OP::ghost();
    m_U0.define(a_map.layout(), ghost);
    m_bc = std::make_shared<BC>(a_map);
}

template <
    template<typename, template<MemType> typename, MemType> class OPType,
    template<MemType> class MAP,
    typename T,
    template<typename, unsigned int, template<MemType> typename, MemType, Centering> class BCType,
    MemType MEM>
void MBLevelRK4<OPType, MAP, T, BCType, MEM>::advance(
    LevelStateData& a_state, 
    LevelScalarData& a_dVolrLev, 
    double&         a_dt,
    double          a_time,
    int             a_order)
{
    double sixth = 1, third=1, half = 1;
    sixth/=6; third/=3; half/=2;
    int kstage;
    m_delta.setVal(0.);

    if (a_order == 1)
    {
        a_state.copyTo(m_U0);
        kstage = 0;
        OP::preStageLevel(m_U0, m_k, a_dVolrLev, m_interp, a_time, a_dt, kstage);
        m_delta.increment(m_k, 1.0);
        OP::StageLevel(m_k, m_U0, a_time, a_dt, kstage);
        m_delta.increment(m_k, 1.0);
    }

    if (a_order == 4)
    {
        a_state.copyTo(m_U0);
        kstage = 0;
        OP::preStageLevel(m_U0, m_k, a_dVolrLev, m_interp, a_time, a_dt, kstage);
        m_delta.increment(m_k, 1.0);
        OP::StageLevel(m_k, m_U0, a_time, a_dt, kstage);
        m_delta.increment(m_k, sixth);

        a_state.copyTo(m_U0);
        m_U0.increment(m_k, half);
        kstage = 1;
        OP::preStageLevel(m_U0, m_k, a_dVolrLev, m_interp, a_time+0.5*a_dt, a_dt, kstage);
        OP::StageLevel(m_k, m_U0, a_time+0.5*a_dt, a_dt, kstage);
        m_delta.increment(m_k, third);

        a_state.copyTo(m_U0);
        m_U0.increment(m_k, half);
        kstage = 2;
        OP::preStageLevel(m_U0, m_k, a_dVolrLev, m_interp, a_time+0.5*a_dt, a_dt, kstage);
        OP::StageLevel(m_k, m_U0, a_time+0.5*a_dt, a_dt, kstage);
        m_delta.increment(m_k, third);

        a_state.copyTo(m_U0);
        m_U0.increment(m_k);
        kstage = 3;
        OP::preStageLevel(m_U0, m_k, a_dVolrLev, m_interp, a_time+a_dt, a_dt, kstage);
        OP::StageLevel(m_k, m_U0, a_time+a_dt, a_dt, kstage);
        m_delta.increment(m_k, sixth);
    }

    a_state.increment(m_delta, 1.0);

}

} //end Proto namespace

#endif //end include guard
