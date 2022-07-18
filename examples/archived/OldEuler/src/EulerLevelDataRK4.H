#ifndef _EULER_LEVEL_DATA_RK4_H_
#define _EULER_LEVEL_DATA_RK4_H_

#include "Proto.H"
#include "Proto_LevelBoxData.H"
#include "EulerOp.H"
#include <iostream>
using namespace Proto;

//#define NCOMPS DIM+2

class EulerLevelDataDX;

class EulerLevelDataState
{
public:
    EulerLevelDataState();
    ~EulerLevelDataState();
    EulerLevelDataState(const ProblemDomain& a_probDom,
                        const Point& a_boxSize,
                        const double a_dx,
                        const double a_gamma);
    double m_dx;
    double m_gamma;
    //TODO: Put in MEMTYPE so we can run on the GPU?
    Reduction<double, Abs> m_Rxn;
    LevelBoxData<double,NUMCOMPS> m_U;
    DisjointBoxLayout m_dbl;
    void increment(const EulerLevelDataDX& a_DX);
};

class EulerLevelDataDX
{
public:
    EulerLevelDataDX();
    ~EulerLevelDataDX();
    void init(EulerLevelDataState& a_State);
    void increment(double& a_weight, const EulerLevelDataDX& a_incr);
    void operator*=(const double& a_weight);
    //TODO: Put in MEMTYPE so we can run on the GPU?
    LevelBoxData<double,NUMCOMPS> m_DU;
    DisjointBoxLayout m_dbl;
};

class EulerLevelDataRK4Op
{
public:
    EulerLevelDataRK4Op();
    ~EulerLevelDataRK4Op();
    void operator()(EulerLevelDataDX& a_DX,
                    double a_time,
                    double a_dt,
                    EulerLevelDataState& a_State);
};
#endif
