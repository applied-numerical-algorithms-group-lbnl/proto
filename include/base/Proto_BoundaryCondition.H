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

    /// Abstract Boundary Condition Interface 
    template<class T, unsigned int C=1, MemType MEM=HOST, unsigned int D=1, unsigned int E=1>
    class BoundaryCondition {
        public:
        typedef BoxData<T,C,MEM,D,E> StateData;

        /// Query Boundary Box
        static inline Box Boundary(const LevelIndex& a_index, Face a_face);
        
        /// Compute Dirichlet Stencil
        static inline BCStencil<T> DirichletStencil(T a_dx, Face a_face);
        
        /// Apply Constant Dirichlet BC to all values
        static inline void Dirichlet(
                BoxData<T,C,MEM,D,E>& a_flux,
                const BoxData<T,C,MEM,D,E>& a_state,
                Box a_domainBox,
                T a_value,
                T a_dx,
                Face a_face);
        
        /// Apply Variable Dirichlet BC
        static inline void Dirichlet(
                BoxData<T,C,MEM,D,E>& a_flux,
                const BoxData<T,C,MEM,D,E>& a_state,
                Box a_domainBox,
                const BoxData<T,C,MEM,D,E>& a_values,
                T a_dx,
                Face a_face);

        /// Base Constructor
        BoundaryCondition(DisjointBoxLayout a_layout, Face a_face)
            : m_layout(a_layout), m_face(a_face) {}

        /// Apply BC Interface
        inline virtual void apply(
                StateData& a_flux,
                const StateData& a_state,
                const LevelIndex& a_index,
                Array<T,DIM> a_dx) const = 0;

        /// Query Face
        inline Face face() const { return m_face; }
        
        /// Query Layout
        inline const DisjointBoxLayout& layout() { return m_layout; }
        protected:
        DisjointBoxLayout m_layout;
        Face m_face;
    };

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

    #include "implem/Proto_BoundaryConditionImplem.H"
} // end namespace Proto
#endif //end include guard
