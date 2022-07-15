#pragma once
#ifndef _PROTO_LEVEL_RK4_H_
#define _PROTO_LEVEL_RK4_H_

#include "ProtoAMR.H"

int LEVELRK4_STAGE;
int LEVELRK4_ITER;

namespace Proto {

///  Generic Explicit RK4 Algorithm
/**  Given y'=f(y,t), it computes y_{t+1} from y_t using the following equations:
     k1=dt*f(y_t,t)
     k2=dt*f(y_t+0.5*k1,t+0.5*dt)
     k3=dt*f(y_t+0.5*k2,t+0.5*dt)
     k4=dt*f(y_t+k3,t+dt)
     y_{t+1}=y_t+k1/6+k2/3+k3/3+k4/6

*/

template <template<typename, MemType> class OPType,
typename T, MemType MEM = MEMTYPE_DEFAULT>
class LevelRK4
{
public:
    typedef OPType<T, MEM> OP; 
    typedef LevelBoxData<T, OP::numState(), MEM, PR_CELL> LevelStateData;
    
    LevelRK4(DisjointBoxLayout& a_layout, T a_dx);
    //Compute y_{t+1} at time a_time and with y_t=a_state, and place it in a_state.
    void advance(LevelStateData& a_state, double& a_dt, double& a_time = 0);
protected:
    LevelStateData m_k;
    LevelStateData m_delta;
    LevelOp<OPType, T, MEM> m_f;
    LevelStateData m_U0;
};

template <template<typename, MemType> class OPType, typename T, MemType MEM>
LevelRK4<OPType, T, MEM>::LevelRK4(
    DisjointBoxLayout& a_layout, T a_dx)
{
    m_f.define(a_dx);
    m_delta.define(a_layout, Point::Zeros());
    m_k.define(a_layout, Point::Zeros());
    m_U0.define(a_layout, OP::ghost());
}

template <template<typename, MemType> class OPType, typename T, MemType MEM>
void LevelRK4<OPType, T, MEM>::advance(
    LevelStateData& a_state, 
    double&         a_dt,
    double&         a_time)
{
    HDF5Handler h5;
    double sixth = 1, third=1, half = 1;
    sixth/=6; third/=3; half/=2;

    m_delta.setVal(0);

    m_f.setTime(a_time);
    a_state.exchange();
    m_f(m_k, a_state, a_dt);    // compute k1*dt
    m_delta.increment(m_k, sixth);
    
    a_state.copyTo(m_U0);
    m_U0.increment(m_k, half);
    m_U0.exchange();
    m_f.setTime(a_time+half*a_dt);
    m_f(m_k, m_U0, a_dt);  // compute k2*dt
    m_delta.increment(m_k, third);
    
    a_state.copyTo(m_U0);
    m_U0.increment(m_k, half);
    m_U0.exchange();
    m_f(m_k, a_state, a_dt);  // compute k3*dt
    m_delta.increment(m_k, third);
    
    m_f.setTime(a_time+a_dt);
    a_state.copyTo(m_U0);
    m_U0.increment(m_k);
    m_U0.exchange();
    m_f(m_k, a_state, a_dt); // compute k4*dt
    m_delta.increment(m_k, sixth);
    
    a_state.increment(m_delta);
}

} //end Proto namespace

#endif //end include guard
