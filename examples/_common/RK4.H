#ifndef _PROTO_RK4_H_
#define _PROTO_RK4_H_
#include "Proto.H"

namespace Proto {

///  Generic Explicit RK4 Algorithm
/**  Given y'=f(y,t), it computes y_{t+1} from y_t using the following equations:
     k1=0.5*dt*f(y_t,t)
     k2=0.5*dt*f(y_t+k1,t+0.5*dt)
     k3=dt*f(y_t+k2,t+0.5*dt)
     k4=dt*f(y_t+k3,t+dt)
     y_{t+1}=y_t+k1/6+k2/3+k3/3+k4/6

     This class is templated on three other classes X, dX, and F.
     These classes have the following requirements:

     //Container for y as well as any parameters that define the application.
     //See EulerState in examples/Euler/src/EulerRK4.H for an example.
     //Note that the Euler parameter gamma and the grid spacing dx are stored in this class.
     class X
     {
       // Parameters and data
     };

     //Helper class to increment X at each intermediate step. It stores the current dX.
     //See EulerDX in examples/Euler/src/EulerRK4.H for an example.
     class dX
     {
       public:
       //Initialize the underlying dX data to a_State
       void init(X& a_State);

       //Performs the operation m_dX+=a_weight*a_incr, where m_dX is the underlying data
       void increment(double& a_weight,const dX& a_incr);

       //Performs the operation m_dX*=a_weight
       void operator*=(const double& a_weight);
     };

     //Helper class to compute f(y,t). This class doesn't contain any state.
     //See EulerOp in examples/Euler/src/EulerRK4.H for an example.
     class F
     {
     public:
     //Computes a_dt*f(a_State+a_DX,a_time) and puts the result in a_DX.
     //NOTE it returns f scaled by dt
     void operator()(dX& a_DX, double a_time, double a_dt, X& a_State);
     };

     WARNINGS:
     1) There is no way to pass parameters or additional information to RK4. All of this information should
     be placed in the state class X.
     2) F computes dt*f(y,t), so make sure to scale by dt.
*/
template <class X, class F, class dX>
class RK4
{
public:
    //Compute y_{t+1} at time a_time and with y_t=a_state, and place it in a_state.
    void advance(double& a_time, double& a_dt, X& a_state);
protected:
    dX m_k;
    dX m_delta;
    F m_f;
};

template <class X, class F, class dX>
void RK4<X, F, dX>::advance(double& a_time, double& a_dt, X& a_state)
{
    double sixth = 1, third=1, half = 1;
    sixth/=6; third/=3; half/=2;

    m_delta.init(a_state);
    m_k.init(a_state);                  // init must allocate stroage, and initialize it to zero.
    m_f(m_k, a_time, a_dt, a_state);    // compute k1
    m_delta.increment(sixth, m_k);
    m_k*=half;
    m_f(m_k, a_time+half*a_dt, a_dt, a_state);  // compute k2
    m_delta.increment(third, m_k);
    m_k*=half;
    m_f(m_k, a_time+half*a_dt, a_dt, a_state);  // conpute k3
    m_delta.increment(third, m_k);
    m_f(m_k, a_time+a_dt, a_dt, a_state); // compute k4
    m_delta.increment(sixth, m_k);
    a_state.increment(m_delta);
}

} //end Proto namespace

#endif //end include guard
