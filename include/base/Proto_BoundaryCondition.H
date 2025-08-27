#pragma once
#ifndef __PROTO_BOUNDARY_CONDITION__
#define __PROTO_BOUNDARY_CONDITION__

#include "Proto.H"
namespace Proto {

    /// @brief A variation of Stencil used for applying boundary conditions. 
    /// Consider the diagram below. BCStencil is designed to extrapolate the
    /// value of a ghost cell (A) using interior data (C) in addition to an edge centered
    /// value on the boundary (B). Note that a different MBStencil is needed to fill
    /// ghost cells at different distances from the boundary. The updated value (A)
    /// may have different units than the source values (B, C, D, ...); for example, A might
    /// be a flux while the source values have units of the associated state.
    /// The diagram below corresponds to the case where a ghost cell in the second layer of 
    /// a low-face is being filled. 
    ///
    ///  ghost <-|-> interior
    ///  +---+---+---+---+- ...
    ///  | A |   B C | D |  ...
    ///  +---+---+---+---+  ...
    ///
    /// As an equation:
    /// A = B*coef + stencil(C, D,...)
    /// 
    /// This structure is used primarily in the implementations to follow and is
    /// not recommended for public use.
    template<class T>
    struct BCStencil {
        Stencil<T> stencil; 
        T coef;
    };

    /// @brief Library of utility functions for assisiting in the application of boundary conditions
    class BoundaryCondition {
        public:

        /// @brief Returns the Box corresponding to the boundary of domainBox on face. 
        /// For low-faces along the axis dir, the output is domainBox.edge(dir).
        /// For high-faces, the output is domainBox.adjacent(dir). That is to say,
        /// Low boundaries are contained within domainBox and high boundaries are just
        /// outside of domainBox. This is the convention used in all of Proto. 
        /// @param domainBox - The Box corresponding to the interior of a domain
        /// @param face - The face of domainBox whose boundary Box is to be computed
        static inline Box Boundary(Box domainBox, Face face);
        
        /// @brief Generates the BCStencil object needed to apply a dirichlet boundary
        /// condition to a flux in direction dir on a specified face
        /// @param dx - Grid spacing
        /// @param face - The face corresponding to the boundary condition
        /// @param fluxDir - Coordinate of the flux being updated. Can be normal or tangential
        template<class T>
        static inline BCStencil<T> DirichletStencil(T dx, Face face, int fluxCoord);
        
        /// @brief Applies a constant dirichlet boundary condition on a specified face by updating
        /// the fluxes in all directions
        /// @param fluxes - flux data which will be updated on the appropriate boundary to apply the BC
        /// @param state - state data used to apply the BC
        /// @param dirichletBCValue - constant value of the state on the boundary. Note that the same constant
        /// Will be applied to every component of each flux. Alias as needed to apply BCs componentwise.
        /// @param face - the face on which the boundary condition will be applied
        /// @param interiorBox - defines which elements of state and fluxes are interior to the domain.
        /// This can be larger than the Boxes defining the state / fluxes. In practice, fluxes and state
        /// are patch scoped values and interiorBox is the Box associated with a ProblemDomain.
        /// @param dx - anisotropic grid spacing
        template<class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
        static inline void Dirichlet(
                Array<BoxData<T,C,MEM,D,E>,DIM>& fluxes,
                const BoxData<T,C,MEM,D,E>& state,
                T dirichletBCValue,
                Face face,
                Box interiorBox,
                Array<T,DIM> dx);
        
        /// @brief Applies a variable dirichlet boundary condition on a specified face by updating
        /// the fluxes in all directions
        /// @param fluxes - flux data which will be updated on the appropriate boundary to apply the BC
        /// @param state - state data used to apply the BC
        /// @param dirichletBCValues - values of the state on the boundary
        /// @param face - the face on which the boundary condition will be applied
        /// @param interiorBox - defines which elements of state and fluxes are interior to the domain.
        /// This can be larger than the Boxes defining the state / fluxes. In practice, fluxes and state
        /// are patch scoped values and interiorBox is the Box associated with a ProblemDomain.
        /// @param dx - anisotropic grid spacing
        template<class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
        static inline void Dirichlet(
                Array<BoxData<T,C,MEM,D,E>,DIM>& fluxes,
                const BoxData<T,C,MEM,D,E>&      state,
                const BoxData<T,C,MEM,D,E>&      dirichletBCvalues,
                Face face,
                Box interiorBox,
                Array<T,DIM> dx);

        #ifdef PR_LAPACK

        /// @brief Creates a BCStencil which is used to extrapolate ghost cell values a specified 
        /// distance from a specified face. The number of source data points can also be specified.
        /// @tparam T 
        /// @param dist - distance of the extrapolated ghost cell from the boundary.
        /// dist = 1 is a cell directly adjacent to the boundary
        /// @param face - the face from which extrapolation is computed
        /// @param npoints - the number of interior data points used to compute the extrapolation.
        template<class T>
        static inline BCStencil<T> DirichletExtrapStencil(int dist, Face face, int npoints);
        
        /// @brief Extrapolates state data into a ghost region in accordance to a specified constant
        /// dirichlet boundary condition. This function should only be used for data which is very smoothly
        /// varying near the boundary in question.  
        /// @param state - data updated by and used as a source for the extrapolation. Elements inside of 
        /// interiorBox are inputs to the extrapolation and elements outside of interiorBox are assumed to be
        /// in the ghost region and will be computed by the extrapolation.
        /// @param dirichletBCValue - constant value of the dirichlet boundary condition
        /// @param face - face on which the dirichlet boundary condition is applied
        /// @param interiorBox - Box which defines which elements of state are in the interior or ghost region
        template<class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
        static inline void DirichletFillGhost(
            BoxData<T,C,MEM,D,E>& state,
            T dirichletBCValue,
            Face face,
            Box interiorBox);

        /// @brief Extrapolates ghost region data in boundary regions of codimension 2 or greater (e.g. hyper-corners).
        /// This function assumes that ghost data in boundaries of codimension 1 have already been computed using e.g. 
        /// DirichletFillGhost or some other boundary filling strategy.
        /// @param state - data updated by and used as a source for the extrapolation. Elements which are in boundary regions
        /// - as defined by interiorBox - of codimension 2 or greater are updated while elements in boundary regions of codimension
        /// 1 are used as source data. Interior data are not used or updated.
        /// @param interiorBox - used to define which elements of state are interior cells, boundary cells, etc.
        template<class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
        static inline void ExtrapolateCorners(
            BoxData<T,C,MEM,D,E>& state,
            Box interiorBox);
        #endif

    };

// The following code is part of an old policy-based implementation wherein boundary conditions were implemented as first class
// objects. This may still be a good approach, but I suspect it's more flexible and user-friendly to require users to manually apply
// their own boundary conditions and supply helpful functions like those above to make this process as concise and expressive as possible.
// - Chris Bozhart - 8/21/2025
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
