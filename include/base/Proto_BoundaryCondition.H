#pragma once
#ifndef __PROTO_BOUNDARY_CONDITION__
#define __PROTO_BOUNDARY_CONDITION__

#include "Proto.H"
namespace Proto {

    template<class T>
    struct BCStencil {
        Stencil<T> stencil;
        T coef;
    };

    class BoundaryCondition {
        public:

        /// Query Boundary Box
        static inline Box Boundary(Box a_domainBox, Face a_face);
        
        /// Compute Dirichlet Stencil
        template<class T>
        static inline BCStencil<T> DirichletStencil(T a_dx, Face a_face, int a_dir);
        
        /// Apply Constant Dirichlet BC to all values
        template<class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
        static inline void Dirichlet(
                Array<BoxData<T,C,MEM,D,E>,DIM>& a_fluxes,
                const BoxData<T,C,MEM,D,E>&      a_state,
                T a_value,
                Face a_face,
                Box a_domainBox,
                Array<T,DIM> a_dx);
        
        /// Apply Variable Dirichlet BC
        template<class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
        static inline void Dirichlet(
                Array<BoxData<T,C,MEM,D,E>,DIM>& a_flux,
                const BoxData<T,C,MEM,D,E>&      a_state,
                const BoxData<T,C,MEM,D,E>& a_values,
                Face a_face,
                Box a_domainBox,
                Array<T,DIM> a_dx);

        #ifdef PR_OPS
        template<class T>
        static inline BCStencil<T> DirichletExtrapStencil(int dist, Face face, int npoints);
        
        template<class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
        static inline void DirichletFillGhost(
            BoxData<T,C,MEM,D,E>& state,
            T value,
            Face face,
            Box domainBox);
        #endif

    };
#if 0
    template<class T, unsigned int C=1, MemType MEM=HOST, unsigned int D=1, unsigned int E=1>
    class BC_ConstDirichlet : public BoundaryCondition<T,C,MEM,D,E> {
        public: 
        typedef BoxData<T,C,MEM,D,E> StateData;

        BC_ConstDirichlet(Face a_face, T a_value, T a_dx)
            : BoundaryCondition<T,C,MEM,D,E>(a_face), m_value(a_value), m_dx(a_dx) {}

        virtual void apply(
                StateData& a_flux,
                const StateData& a_state,
                Box a_domainBox) const
        {
            this->Dirichlet(a_flux, a_state, a_domainBox, m_value, m_dx, this->face());
        }

        T& data() { return m_value; }

        private:
        T m_value;
        const T m_dx;
    };
    
    template<class T, unsigned int C=1, MemType MEM=HOST, unsigned int D=1, unsigned int E=1>
    class BC_VariableDirichlet : public BoundaryCondition<T,C,MEM,D,E> {
        public: 
        typedef BoxData<T,C,MEM,D,E> StateData;

        BC_VariableDirichlet(Face a_face, BoxData<T,C,MEM,D,E>& a_values, T a_dx)
            : m_values(a_values), m_dx(a_dx) { this->m_face = a_face; }

        virtual void apply(
                StateData& a_flux,
                const StateData& a_state,
                Box a_domainBox) const
        {
            this->Dirichlet(a_flux, a_state, a_domainBox, m_values, m_dx, this->face());
        }

        BoxData<T,C,MEM,D,E>& data() { return m_values; }

        private:
        BoxData<T,C,MEM,D,E>& m_values;
        const T m_dx;
    };
    #endif
    #include "implem/Proto_BoundaryConditionImplem.H"
} // end namespace Proto
#endif //end include guard
