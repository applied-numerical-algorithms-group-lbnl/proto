#ifndef _MHD_LEVEL_DATA_RK4_H_
#define _MHD_LEVEL_DATA_RK4_H_

#include "Proto.H"
#include "Proto_LevelBoxData.H"
#include "MHDOp.H"
#include <iostream>
using namespace Proto;

#define NCOMPS 2*DIM+2

class MHDLevelDataDX;

class MHDLevelDataState
{
public:
    MHDLevelDataState();
    ~MHDLevelDataState();
    MHDLevelDataState(const ProblemDomain& a_probDom,
                        const Point& a_boxSize,
                        const double a_dx,
                        const double a_gamma);
    double m_dx;
    double m_gamma;
    //TODO: Put in MEMTYPE so we can run on the GPU?
    LevelBoxData<double,NUMCOMPS> m_U;
    DisjointBoxLayout m_dbl;
    void increment(const MHDLevelDataDX& a_DX);
};

class MHDLevelDataDX
{
public:
    MHDLevelDataDX();
    ~MHDLevelDataDX();
    void init(MHDLevelDataState& a_State);
    void increment(double& a_weight, const MHDLevelDataDX& a_incr);
    void operator*=(const double& a_weight);
    //TODO: Put in MEMTYPE so we can run on the GPU?
    LevelBoxData<double,NUMCOMPS> m_DU;
    DisjointBoxLayout m_dbl;
};

class MHDLevelDataRK4Op
{
public:
    MHDLevelDataRK4Op();
    ~MHDLevelDataRK4Op();
    void operator()(MHDLevelDataDX& a_DX,
                    double a_time,
                    double a_dt,
                    MHDLevelDataState& a_State);
};

class MHDLevelDataEulerOp
{
public:
    MHDLevelDataEulerOp();
    ~MHDLevelDataEulerOp();
    void operator()(MHDLevelDataDX& a_DX,
                    double a_time,
                    double a_dt,
                    MHDLevelDataState& a_State);
};

class MHDLevelDataViscosityOp
{
public:
    MHDLevelDataViscosityOp();
    ~MHDLevelDataViscosityOp();
    void operator()(MHDLevelDataDX& a_DX,
                    double a_time,
                    double a_dt,
                    MHDLevelDataState& a_State);
};
#endif
