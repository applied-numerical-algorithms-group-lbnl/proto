#pragma once

#include "IMEXIntegrator.H"

/// Midpoint Method IMEX ARK Integrator
/**
    Simple IMEX ARK integrator using the Midpoint Method.
*/

template<class IMEXOP, class IMEXDATA>
class Midpoint_IMEXIntegrator : public IMEXIntegrator<IMEXOP, IMEXDATA, 2,1>
{
    void define() {
        this-> m_aIDiag = 1.0;

        // Define 1D arrays
        this->m_c = {0.0,0.5};
        this->m_b = {0.0,1.0};

        // Define MD arrays
        this->m_bstar = { {0.0,1.0} };
        this->m_aE = {{0.0, 0.0},
            {0.5, 0.0}};
        this->m_aI = {{0.0, 0.0},
            {0.0, 0.5}};

        this->m_defined = true;
    }
};
