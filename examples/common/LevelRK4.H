#pragma once
#ifndef _PROTO_LEVEL_RK4_H_
#define _PROTO_LEVEL_RK4_H_

#include "ProtoAMR.H"

namespace Proto {

///  Generic Explicit RK4 Algorithm
/**  Given y'=f(y,t), it computes y_{t+1} from y_t using the following equations:
     k1=0.5*dt*f(y_t,t)
     k2=0.5*dt*f(y_t+k1,t+0.5*dt)
     k3=dt*f(y_t+k2,t+0.5*dt)
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
    LevelStateData m_U0;
    LevelStateData m_k;
    LevelStateData m_delta;
    LevelOp<OPType, T, MEM> m_f;
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
    m_f.setFluxScale(a_dt);
    m_f.setTime(a_time);
    a_state.copyTo(m_U0);
    
    m_f(m_k, m_U0);    // compute k1
    /*
    h5.writeLevel(m_f.dx()[0], m_k, "K1");
    m_delta.increment(m_k, sixth);
    m_k*=half;
    m_f.setTime(a_time+half*a_dt);
    m_f(m_k, m_U0);  // compute k2
    h5.writeLevel(m_f.dx()[0], m_k, "K2");
    m_delta.increment(m_k, third);
    m_k*=half;
    m_f(m_k, m_U0);  // conpute k3
    h5.writeLevel(m_f.dx()[0], m_k, "K3");
    m_delta.increment(m_k, third);
    m_f.setTime(a_time+a_dt);
    m_f(m_k, m_U0); // compute k4
    h5.writeLevel(m_f.dx()[0], m_k, "K4");
    m_delta.increment(m_k, sixth);
    a_state.increment(m_delta);
    */
}

} //end Proto namespace

#endif //end include guard
