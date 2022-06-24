#ifndef _LEVELMULTIGRID_H_
#define _LEVELMULTIGRID_H_
#include <cmath>
#include <cstring>
#include <memory>
#include "Proto.H"
#define MAXBOXSIZE 64
using namespace std;
using namespace Proto;

/// Multigrid: solve a FV discretization of Poisson's equation on a rectangle.
class LevelMultigrid
{
    public:
        LevelMultigrid(
                const DisjointBoxLayout& a_layout,
                double a_dx,
                int a_level);

        void define(
                const DisjointBoxLayout& a_layout,
                double a_dx,
                int a_level);

        void coarseResidual(
                LevelBoxData<double >& a_resc,
                LevelBoxData<double >& a_phi,
                LevelBoxData<double >&a_rhs);
        double resnorm(
                LevelBoxData<double >& a_phi,
                LevelBoxData<double >& a_rhs);

        void pointRelax(
                LevelBoxData<double >& a_phi,
                LevelBoxData<double >& a_rhs,
                int a_numIter);

        void avgDown(
                LevelBoxData<double >& a_resc,
                LevelBoxData<double >& a_res);

        void fineInterp(
                LevelBoxData<double >& a_phiFine,
                LevelBoxData<double >& a_deltaCoarse);

        void vCycle(
                LevelBoxData<double >& a_phi,
                LevelBoxData<double >& a_rhs );
    private:
        LevelBoxData<double > m_localCoarse;
        LevelBoxData<double > m_resc;
        LevelBoxData<double > m_delta;
        Box m_box;
        DisjointBoxLayout m_layout;
        int m_domainSize;
        shared_ptr<LevelMultigrid > m_coarsePtr;
        double m_dx;
        double m_lambda;
        int m_level;
        int m_preRelax = 2*DIM;
        int m_postRelax = 2*DIM;
        int m_bottomRelax = 10;
        InterpStencil<double> m_fineInterp = InterpStencil<double>::PiecewiseConstant(Point::Ones()*2); 
        Reduction<double, Abs> m_rxn;
};
#endif
