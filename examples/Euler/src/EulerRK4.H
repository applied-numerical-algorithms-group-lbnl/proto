#ifndef _EULERRK4_H_
#define _EULERRK4_H_

#include "EulerOp.H"
#include "Proto.H"

using namespace Proto;

class EulerDX;
class EulerState
{
public:
    EulerState(){};
    ~EulerState(){};
    EulerState(const Bx& a_box);
    BoxData<double,NUMCOMPS> m_U;
    Bx m_dbx0;
    double m_velSave;
    void increment(const EulerDX& a_DX);
};

class EulerDX
{
public:
    EulerDX(){};
    ~EulerDX(){};
    BoxData<double,NUMCOMPS> m_DU;
    Bx m_box;
    void increment(double& a_weight,const EulerDX& a_incr);
    void init(const EulerState& a_State);
    void operator*=(const double& a_weight);
};

class EulerRK4Op
{
public:
    EulerRK4Op(){};
    ~EulerRK4Op(){};
    static int s_count;
    EulerOp m_eop;
    void operator()(
                  EulerDX& a_DX,
                  double a_time,
                  double a_dt,
                  EulerState& a_State);
};
#endif
