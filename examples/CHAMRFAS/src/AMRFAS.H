#pragma once

#include "Multigrid.H"
#include "Proto.H"
#include "NeighborIterator.H"
#include "AMRUtils.H"
#include "LevelFluxRegister.H"

//using namespace Proto;

template<class AMR_OP, class DATA = BaseFab<Real>>
class AMRFAS
{
    public:
        AMRFAS(){m_defined = false;}

        AMRFAS(const std::vector<DisjointBoxLayout>& a_grids,
                double a_dx,
                int    a_level,
                int    a_bottomDepth)
        {
            define(a_grids, a_dx, a_level, a_bottomDepth);
        }
        inline void define(const std::vector<DisjointBoxLayout>& a_grids,
                double a_dx,
                int    a_level,
                int    a_bottomDepth);

        inline void vcycle(std::vector<std::shared_ptr<LevelData<DATA>>> a_phi,
                const std::vector<std::shared_ptr<LevelData<DATA>>> a_rho,
                const std::vector<std::shared_ptr<LevelData<DATA>>> a_res);

        inline void vcycle(std::vector<std::shared_ptr<LevelData<DATA>>> a_phi,
                const std::vector<std::shared_ptr<LevelData<DATA>>> a_rho,
                const std::vector<std::shared_ptr<LevelData<DATA>>> a_res,
                const LevelData<DATA>& a_R);

        inline void operator()(std::vector<std::shared_ptr<LevelData<DATA>>> a_rhs,
                std::vector<std::shared_ptr<LevelData<DATA>>> a_phi);
        template<typename... Args>
        inline void write(std::vector<std::shared_ptr<LevelData<FArrayBox>>> a_data,
                const char* a_fname,
                Args... a_params);

    private:
        int                        m_level;

        Multigrid<AMR_OP, DATA>    m_mg;
        AMR_OP                     m_op; 
        LevelData<DATA>            m_phiC0;
        LevelData<DATA>            m_RC;

        std::shared_ptr<AMRFAS<AMR_OP, DATA>> m_coarser;
        LevelFluxRegister          m_reflux;
        bool m_defined;

        // implementation of course-fine boundary
};

#include "AMRFASImplem.H"
