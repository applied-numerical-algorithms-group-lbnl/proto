#ifndef _MHD_LEVEL_DATA_RK4_H_
#define _MHD_LEVEL_DATA_RK4_H_

#include "Proto.H"
#include "Proto_LevelBoxData.H"

constexpr unsigned int NUMCOMPS=2*DIM+2;

#define NGHOST 6

#define PI 3.1415926535897932384626433832795


using namespace Proto;

typedef Var<double,DIM> V;
typedef Var<double,NUMCOMPS> State;

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
	                  const double a_dy,
	                  const double a_dz,
	                  const double a_gamma);
	double m_dx;
	double m_dy;
	double m_dz;
	double m_gamma;
//TODO: Put in MEMTYPE so we can run on the GPU?
	LevelBoxData<double,NUMCOMPS> m_U;
	DisjointBoxLayout m_dbl;
	ProblemDomain m_probDom;
	void increment(const MHDLevelDataDX& a_DX);
};
#endif
