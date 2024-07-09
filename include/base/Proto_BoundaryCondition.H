#pragma once
#ifndef __PROTO_BOUNDARY_CONDITION__
#define __PROTO_BOUNDARY_CONDITION__

#include "Proto.H"
namespace Proto {
    
    /// Abstract Boundary Condition Interface 
    template<class T>
    class BoundaryCondition {
        public:
        
        template<unsigned int C, MemType MEM, unsigned int D, unsigned int E>
        virtual void apply(
                BoxData<T,C,MEM,D,E>& a_flux,
                const BoxData<T,C,MEM,D,E> a_state,
                Box a_domainBox) const = 0;

        Face face() const { return m_face; }

        Box boundary(Box a_domainBox) const;
        protected:
        
        const Face m_face;
    };

    template<class T>
    class ConstDirichletBC : public BoundaryCondition<T> {
        public: 

        ConstDirichletBC(T a_value, Face a_face, T a_dx);

        template<unsigned int C, MemType MEM, unsigned int D, unsigned int E>
        void apply(
                BoxData<T,C,MEM,D,E>& a_flux,
                const BoxData<T,C,MEM,D,E> a_state,
                Box a_domainBox) const;
        
        private: 

        const T m_value;
        Stencil<T> m_stencil;
        T m_offset;

    };
    
    template<class T>
    class ConstNeumannBC : public BoundaryCondition {
        public: 

        ConstNeumannBC(Face a_face, T a_value);

        template<unsigned int C, MemType MEM, unsigned int D, unsigned int E>
        void apply(BoxData<T,C,MEM,D,E> a_flux, Box a_domainBox) const;
        
        private: 

        const T m_value;
        Stencil<T> m_stencil;

    };

    #include "implem/Proto_BoundaryConditionImplem.H"
} // end namespace Proto
#endif //end include guard
