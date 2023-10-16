#pragma once

#include "IMEXIntegrator.H"

/// 2nd Order Forward - Backward Euler IMEX ARK Method
/**
    Simple IMEX ARK method using forward and backward Euler.
*/

template<class IMEXOP, class IMEXDATA>
class Euler_IMEXIntegrator : public IMEXIntegrator<IMEXOP, IMEXDATA, 2,1>
{
    void define() {
        this-> m_aIDiag = 1.0;

        // Define 1D arrays
        this->m_c = {0.0,1.0};
        this->m_b = {0.0,1.0};

        // Define MD arrays
        this->m_bstar = { {0.0,1.0} };
        this->m_aE = {{0.0, 0.0},
            {1.0, 0.0}};
        this->m_aI = {{0.0, 0.0},
            {0.0, 1.0}};

        this->m_defined = true;
    }
};
