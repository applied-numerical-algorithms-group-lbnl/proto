#pragma once

#include "AMRFAS.H"
#include "IMEXOp.H"

/// AMRFAS IMEX Operator
/**
    Operator class compatable with the IMEX framework which uses the AMRFAS framework for its solver.

    \tparam IMOP    The implicit operator physics; this is a class that extends BaseOp
    \tparam EXOP    The explicit operator physics; this is a class that extends BaseOp
    \tparam NCOMP   Number of components
*/
template<class IMOP, class EXOP, int NCOMP=1>
class AMRFAS_IMEXOp : IMEXOp<AMRData<NCOMP>>
{
    public:
    
    void define(AMRLayout& a_layout, int a_bottomDepth)
    {
        m_exop.define(a_layout, this->m_dt, a_bottomDepth);
        m_imop.define(a_layout, this->m_dt, a_bottomDepth);
    }

    void explicitOp(AMRData<NCOMP>& a_result, AMRData<NCOMP>& a_state, Real a_time, int a_stage)
    {
        for (int ii = 0; ii < m_numIterE; ii++)
        {
            m_exop.vcycle(a_state, a_result)
        }
    }
    void implicitOp(AMRData<NCOMP>& a_result, AMRData<NCOMP>& a_state, Real a_time, int a_stage)
    {
        for (int ii = 0; ii < m_numIterI; ii++)
        {
            m_imop.vcycle(a_state, a_result)
        }

    }
    void solve(AMRData<NCOMP>& a_soln, AMRData<NCOMP>& a_rhs, Real a_time, int a_stage) 
    {
        // not yet implemented
    };

    private:
    
    int m_numIterE = 10; // replace this with a convergence check
    int m_numIterI = 10; // replace this with a convergence check

    AMRFAS<EXOP> m_exop;
    AMRFAS<IMOP> m_imop;
};
