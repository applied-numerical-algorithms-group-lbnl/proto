#include "MHDLevelDataRK4.H"


MHDLevelDataState::MHDLevelDataState()
{
}

MHDLevelDataState::~MHDLevelDataState()
{
}

MHDLevelDataState::MHDLevelDataState(const ProblemDomain& a_probDom,
                                     const Point& a_boxSize,
                                     const double a_dx,
                                     const double a_dy,
                                     const double a_dz,
                                     const double a_gamma) :
    m_dx(a_dx),
    m_dy(a_dy),
    m_dz(a_dz),
	m_gamma(a_gamma),
	m_dbl(a_probDom, a_boxSize),
	m_probDom(a_probDom)
	//m_U(m_dbl,Point::Zero())
{
	m_U.define(m_dbl,Point::Zero());
	
}
