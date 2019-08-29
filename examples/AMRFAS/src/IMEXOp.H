#pragma once

#include "Proto.H"

/// IMEX Operator Interface
/**
    This abstract class defines the interface a user must implement to be compatable with the IMEXIntegrator Framework
    TODO: Quit being a lazy clown and add some math to the documentation of this class.

    \tparam IMEXDATA    A data holder class that implements the interface in IMEXData (but may not necessarily inherit from it)
*/
template<class IMEXDATA>
class IMEXOp {

    public:

    /// Default Constructor    
    IMEXOp()
    {
        std::vector<double> dx;
        dx.resize(DIM, 1);
        init(1,1,dx);
    }

    /// Constructor
    IMEXOp(double a_dt, double a_dtScale, std::vector<double>& a_dx)
    {
        init(a_dt, a_dtScale, a_dx);
    }

    /// Lazy Constructor
    void init(double a_dt, double a_dtScale, std::vector<double>& a_dx)
    {
        m_dt = a_dt;
        m_dtScale = a_dtScale;
        for (int ii = 0; ii < DIM; ii++)
        {
            m_dx[ii] = a_dx[ii];
        }
        define();
    }

    /// Set Timestep
    void setTimestep(double a_dt){ m_dt = a_dt; }
    
    /// Set Grid Spacing
    void setGridSpacing(int a_dir, double a_dx){ m_dx[a_dir] = a_dx; } 
    
    /// Set All Grid Spacings
    void setGridSpacing(double a_dx)
    {
        for (int ii = 0; ii < DIM; ii++)
        {
            m_dx[ii] = a_dx;
        }
    }

    /// Pre Process
    /**
        This function is called on the state before the first stage of the RK method.
        By default this function does nothing.
    
        \param a_data   The state being integrated
    */
    virtual void preProcess(IMEXDATA& a_data){}

    /// Explicit Operator
    /**
        Compute the explicit terms of the right hand side.

        \param a_result     Data holder for the output.
        \param a_state      The input state
        \param a_time       Time of evaluation (may be unused if there is no explicit time dependence)
        \param a_stage      Stage number (often unused)
    */
    virtual void explicitOp(IMEXDATA& a_result, IMEXDATA& a_state, double a_time, int a_stage) = 0;
    
    /// Implicit Operator
    /**
        Compute the implicit terms of the right hand side.

        \param a_result     Data holder for the output.
        \param a_state      The input state
        \param a_time       Time of evaluation (may be unused if there is no explicit time dependence)
        \param a_stage      Stage number (often unused)
    */
    virtual void implicitOp(IMEXDATA& a_result, IMEXDATA& a_state, double a_time, int a_stage) = 0;
    
    /// Implicit Operator
    /**
        Performs the implicit solve (see Ascher et.al. 1997 Eqn 34). TODO: transcribe this equation. 

        \param a_result     Data holder for the output.
        \param a_state      The input state
        \param a_time       Time of evaluation (may be unused if there is no explicit time dependence)
        \param a_stage      Stage number (often unused)
    */
    virtual void solve(IMEXDATA& a_soln, IMEXDATA& a_rhs, double a_time, int a_stage) = 0;
    
    /// Post Process
    /**
        This function is called after the last stage of the RK method. By default it does nothing. 
    */
    virtual void postProcess(IMEXDATA& a_data){}

    /// User Defined Construction
    /**
        This is where a user can define necessary Proto::Stencil objects, numerical parameters, etc.
        It is called automatically by init() and the non-trivial constructor so the user doesn't have to explicitly call this. 
        By default, this function does nothing. 
    */
    virtual void define(){}

    protected:
    
    double m_dt;
    double m_dtScale; // scale on dt during the solve
    double m_dx[DIM];
};
