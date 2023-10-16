#pragma once

#include "IMEXOp.H"

/// Generic Implicit - Explicit Additive Runge Kutta Integrator
/**
    Implements the titular integration scheme with the following simplifying assumptions:
    1.) Explicit and Implicit butcher tables share b and c vectors (See Ascher et. al. 1997)
    2.) Method is Diagonally Implicit

    IMEXIntegrator implements the integrator itself, including space for temporaries (k values) and the Butcher Tableau
    Specific integrators (ARK4(3)6L, Forward-Backward Euler, Classic RK4) are implemented as subclasses of IMEXIntegrator.
    Dense output is also implemented (TODO: OK not yet, but it will be!)

    \tparam IMEXOP      A class that extends IMEXOp and implements the interface therein. This is where the physics live.
    \tparam IMEXDATA    A data holder class which implements the interface described in IMEXData.H, but does not HAVE to extend IMEXData.
    \tparam NSTAGES     The number of Runge Kutta stages
    \tparam NDENSE      The number of coefficients for dense output
*/

template<class IMEXOP, class IMEXDATA, int NSTAGES, int NDENSE>
class IMEXIntegrator
{
    public:
    
    /// Default Constructor
    IMEXIntegrator() {
       m_defined = false;
    }

    /// User Specified Define
    /**
        Initializes the Butcher Tableau and other related fields. Defined by the extending class.
        The user does NOT have to call this; it is called automatically by init.
    */
    virtual void define() = 0;

    /// Initialization
    /**
        Initializes all of the temporary data holders using a_state as a template. Also calls define().

        \param a_state          A data holder representing the state being integrated. All data holders in this will be initialized as clones of a_state.
        \param a_dt             Timestep
        \param a_dx             Grid spacings
        \param a_denseOutput    True if the integrator is to compute values for dense output
    */
    virtual void init(IMEXDATA& a_state, double a_dt,
        const std::vector<double>& a_dx, bool a_denseOutput = false);
    /// Advance
    /**
        Advance a_state from time a_time to time a_time + m_dt

        \param a_time   Time before integration
        \param a_state  The state being integrated
    */
    void advance(double a_time, IMEXDATA& a_state);
    
    protected:

    bool m_init = false;
    bool m_defined = false;
    bool m_denseOutput = false;
    bool m_hasDenseOutput = false;
    
    // THIS STUFF IS MANAGED BY INIT(...)
    double m_dt;        ///< Timestep
    double m_dx[DIM];   ///< Grid Spacings
    std::vector<std::shared_ptr<IMEXDATA>> m_phi;   ///< Accumulated state at each stage
    IMEXDATA  m_rhs;                                ///< Temporary RHS
    std::vector<std::shared_ptr<IMEXDATA>> m_denseCoefs; ///< Dense Output Coefficients
    IMEXDATA  m_kE; ///< Temporary for Explicit k stage
    IMEXDATA  m_kI; ///< Temporary for Implicit k stage
    IMEXOP    m_op; ///< IMEX Operator containing the functions for explicit / implicit operators

    
    // ALL THIS STUFF HAS TO BE SET IN DEFINE()
    double m_aIDiag;                          ///<
    std::vector<double> m_c;                  ///< c vector, m_nStages long
    std::vector<double> m_b;                  ///< b vector, m_nStages long
    std::vector<std::vector<double>> m_aE;    ///< Explicit a matrix, m_nStages X m_nStages
    std::vector<std::vector<double>> m_aI;    ///< Implicit a matrix, m_nStages X m_nStages
    std::vector<std::vector<double>> m_bstar; ///< b_star matrix for dense output, m_nDenseCoefs X m_nStages

};

#include "IMEXIntegratorI.H"

