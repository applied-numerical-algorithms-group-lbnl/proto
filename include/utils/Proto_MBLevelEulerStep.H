#pragma once
#ifndef _PROTO_MB_LEVEL_EULERSTEP_H_
#define _PROTO_MB_LEVEL_EULERSTEP_H_

#include "Proto.H"

namespace Proto {

///  Generic Explicit EulerStep Algorithm
/**  Given y'=f(y,t), it computes y_{t+1} from y_t using the following equations:
     k1=dt*f(y_t,t)
     y_{t+1}=y_t+k1
*/

template <
    template<typename, template<MemType> typename, MemType> class OPType,
    template<MemType> class MAP,
    typename T,
    template<typename, unsigned int, template<MemType> typename, MemType, Centering> class BCType = NullBC,
    MemType MEM = MEMTYPE_DEFAULT>
class MBLevelEulerStep
{
public:
    typedef OPType<T, MAP, MEM> OP; 
    typedef BCType<T,OP::numState(), MAP, MEM, PR_CELL> BC;
    typedef MBLevelOp<OPType, MAP, double, BCType, MEM> LOP;
    typedef MBLevelBoxData<T, OP::numState(), MEM, PR_CELL> LevelStateData;
    
    MBLevelEulerStep(
            const MBLevelMap<MAP>&      a_map,
            const MBInterpOp&           a_interp);

    //Compute y_{t+1} at time a_time and with y_t=a_state, and place it in a_state.
    inline void advance(LevelStateData& a_state, double& a_dt, double a_time = 0.);
    inline LOP& op() { return m_f; }
    inline BC& bc() {return *m_bc; }
protected:
    LevelStateData      m_k;
    LevelStateData      m_delta;
    LOP                 m_f;
    std::shared_ptr<BC> m_bc;
    LevelStateData      m_U0;
    const MBInterpOp&   m_interp;
};
template <
    template<typename, template<MemType> typename, MemType> class OPType,
    template<MemType> class MAP,
    typename T,
    template<typename, unsigned int, template<MemType> typename, MemType, Centering> class BCType,
    MemType MEM>
MBLevelEulerStep<OPType, MAP, T, BCType, MEM>::MBLevelEulerStep(
    const MBLevelMap<MAP>&      a_map,
    const MBInterpOp&           a_interp) : m_interp(a_interp)
{
    Array<Point,DIM+1> ghost;
    ghost.fill(Point::Zeros());
    m_f.define(a_map);
    m_delta.define(a_map.layout(), ghost);
    m_k.define(a_map.layout(), ghost);
    ghost.fill(OP::ghost());
    m_U0.define(a_map.layout(), ghost);
    m_bc = std::make_shared<BC>(a_map);
}

template <
    template<typename, template<MemType> typename, MemType> class OPType,
    template<MemType> class MAP,
    typename T,
    template<typename, unsigned int, template<MemType> typename, MemType, Centering> class BCType,
    MemType MEM>
void MBLevelEulerStep<OPType, MAP, T, BCType, MEM>::advance(
    LevelStateData& a_state, 
    double&         a_dt,
    double          a_time)
{
    m_delta.setVal(0);
    m_f.setTime(a_time);
    m_bc->setTime(a_time);
    m_bc->setRKStage(0);
    m_bc->apply(a_state);
    a_state.exchange();
    m_interp.apply(a_state, a_state);
    m_f(m_k, a_state, a_dt);    // compute k1*dt
    m_delta.increment(m_k, 1.0);
    a_state.increment(m_delta);
}

} //end Proto namespace

#endif //end include guard
