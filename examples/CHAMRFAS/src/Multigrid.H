#pragma once

#include "AMRUtils.H"
#include "Proto.H"
#include "LevelData.H"
#include "AMRIO.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"

#include <vector>
#include <iostream>

//using namespace Proto;

template <class OP, class DATA = BaseFab<Real>>
class Multigrid
{
    public:
        typedef LevelData<DATA> LD;

        inline Multigrid(){m_defined = false;}

        inline Multigrid(const DisjointBoxLayout& a_layout,
                double     a_dx,
                int        a_level,
                bool       a_amr) 
        {
            define(a_layout, a_dx, a_level, a_amr);
        };

        inline void define(const DisjointBoxLayout& a_layout,
                double     a_dx,
                int        a_level,
                bool       a_amr);

        // standard MG vcycle with no AMR
        inline void vcycle(LD& a_phi, const LD& a_R);

        // vcycle used with AMR; interpolates boundary conditions
        inline void vcycle(LD& a_phi, const LD& a_phiC, const LD& a_R);

        inline OP& op(){return m_op;}
        void clear();
    private:

        int m_depth;
        int m_level;
        OP  m_op;

        LevelData<DATA> m_phiC0;
        LevelData<DATA> m_phiC;
        LevelData<DATA> m_RC;
        std::shared_ptr<Multigrid<OP,DATA>> m_coarser;
        
        Proto::InterpStencil<Real>  m_amrInterp;
        LevelData<DATA> m_phiCAMR;
        bool m_defined;

};

#include "MultigridImplem.H"
