#pragma once
#ifndef __PROTO_LEVEL_BC__
#define __PROTO_LEVEL_BC__

#include "Proto_LevelBoxData.H"

namespace Proto {

    // [coord, side]
    typedef std::tuple<int, Side::LoHiSide> Face;
    
    public enum BCType
    {
        DIRICHLET_CONST,
        DIRICHLET_ANALYTIC,
        DIRICHLET_EMPIRICAL
        NEUMANN_CONST,
        NEUMANN_ANALYTIC,
        NEUMANN_EMPIRICAL
    }

    /// Abstract class for defining a BC operator
    /**
     * MMB will be handled in a separate MBLevelBCOp 
     *
     * Different BCType will need different constructor signatures
     * */
    template<BCType>
    class LevelBCOp {
        public:


        private:
        Face m_face;
    };
    
    /// Storage for the LevelBCOps on each face 
    class LevelBC {
        public:

        

        private:

        std::map<Face, LevelBCOp> m_ops
    };
    #include "implem/Proto_LevelBCImplem.H"
} // end namespace Proto
#endif //end include guard
