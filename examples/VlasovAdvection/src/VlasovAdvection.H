#ifndef _VLASOV_ADVECTION_H_
#define _VLASOV_ADVECTION_H_

#include "Proto.H"
#include "Proto_LevelBoxData.H"
using namespace Proto;

/**
 * @brief Computes cell-averages of the divergence of the flux N^T*u*phi (put Latex equations here), and places it in a_rhs.
 * @param a_phi Cell-averages of the distribution function
 * @param a_rhs Cell-averages of flux divergences
 * @param a_dx grid spacing in physical space
 * The algorithm is fourth-order accurate, and detailed in the Implementation section of the COGENT Proto 2020 LaTex document.
 * TODO: Put Latex equations here (and don't refer to the document).
 */
void computeRHS(LevelBoxData<double,1>& a_rhs,
                const LevelBoxData<double,1>& a_phi,
                const double a_dx);

#endif
