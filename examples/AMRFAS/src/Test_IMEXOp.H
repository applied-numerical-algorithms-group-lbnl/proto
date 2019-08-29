#pragma once

#include "AMRData.H"
#include "IMEXOp.H"

/// IMEX Framework Test Function
/**
    Simple IMEXOp used for testing
*/
class Test_IMEXOp : public IMEXOp<AMRData<1>>
{
    public:
    
    static constexpr double s_cE = 0;
    static constexpr double s_cI = -0.1;
    static constexpr int NCOMP = 1;
     
    void explicitOp(AMRData<1>& a_result,
        AMRData<1>& a_state,
        double a_time, int a_stage)
    {
        //double coef = s_cE + 1.0/(1.0 + a_time);
        double coef = s_cE;
        a_state.copyTo(a_result);
        a_result.times(coef);
    }
    void implicitOp(AMRData<1>& a_result,
        AMRData<1>& a_state,
        double a_time, int a_stage)
    {
        a_state.copyTo(a_result);
        a_result.times(s_cI);

    }
    void solve(AMRData<1>& a_soln,
        AMRData<1>& a_rhs,
        double a_time, int a_stage)
    {
        double coef = 1.0/(1.0 - m_dtScale*s_cI*m_dt);
        a_rhs.copyTo(a_soln);
        a_soln.times(coef);
    }

    private:
};

