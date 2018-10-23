#pragma once

#include "Multigrid.H"
#include "Proto.H"
#include "NeighborIterator.H"
#include "AMRUtils.H"

using namespace Proto;

template<class AMR_OP, class DATA = BaseFab<Real>>
class AMRFAS
{
public:
    AMRFAS(){m_defined = false;}
    /*
    AMRFAS(const std::vector<DisjointBoxLayout>& a_grids,
           double a_dx,
           int    a_level,
           int    a_bottomDepth)
    {
        //define(a_grids, a_dx, a_level, a_bottomDepth);
    }
    inline void define(const std::vector<DisjointBoxLayout>& a_grids,
                       double a_dx,
                       int    a_level,
                       int    a_bottomDepth);
    
    inline void vcycle(std::vector<LevelData<DATA>>& a_phi,
                       const std::vector<LevelData<DATA>>& a_rho,
                       const LevelData<DATA> a_R);
    
    inline void computeRhs(LevelData<DATA>& a_RC, 
                           const LevelData<DATA>& a_phi,
                           const LevelData<DATA>& a_phiC,
                           const LevelData<DATA>& a_phiCC,
                           const LevelData<DATA>& a_rhoC);
    */
    void buildcf(DisjointBoxLayout& a_fine, DisjointBoxLayout& a_coarse);

private:
    /*
    int                        m_level;
    
    Multigrid<AMR_OP, DATA>    m_mg; 
    AMR_OP&                    m_op; //alias of m_mg.m_op 
    LevelData<DATA>            m_phiC0;
    LevelData<DATA>            m_RC;
     
    std::shared_ptr<AMRFAS<AMR_OP, DATA>> m_coarser;
    LevelFluxRegister          m_reflux;
    */
    bool m_defined;

    // implementation of course-fine boundary
};

#include "AMRFASImplem.H"
