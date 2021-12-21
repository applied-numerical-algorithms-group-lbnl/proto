#pragma once

#include "AMRData.H"
#include "LevelFluxRegister.H"
#include "Proto_WriteBoxData.H"
#include <cmath>

/// Base Operator
/**
    Base class for all operators used in the Proto AMRFAS framework. Contains all necessary utility functions for running Multigrid or AMRMultigrid.
    To create a new operator, the user write a new class that inherits from BaseOp and implements the init() and flux(...) functions specified herein.
    Note that BaseOp - and all classes derived from it - is a level operator; an AMR or Multigrid framework will have a separate BaseOp instance on each level.

    This class implements all of the common computations that are necessary for the Multigrid and AMRFAS algorithms. Hence, this is mostly an internal class in the framework.

    \tparam NCOMP   Number of components in the independent variable
    \tparam NGHOST  Number of layers of ghost cells needed by the operator
*/

template<unsigned int NCOMP, unsigned int NGHOST>
class BaseOp
{
    public:

        typedef Proto::BoxData<Real, NCOMP> patch;
        typedef Proto::Var<Real, NCOMP> var;

        /// Constructor
        /**
            Initializes an undefined BaseOp. Call BaseOp::define before using.
        */
        inline BaseOp() : m_defined(false) {}

        /// Define
        /**
            Weak constructor. Must be called before *this can be used.
            
            \param a_layout     A Chombo DisjointBoxLayout reference. Defines the domain of this level
            \param a_dx         Grid spacing on this level
            \param a_amr        Flag designating whether or not *this is a member of an AMR object (as opposed to a vanilla Multigrid)
        */
        inline void define(const DisjointBoxLayout& a_layout,
                Real a_dx, bool a_amr = true);
        
        /// User Initialization
        /**
            Part of the weak construction that is defined by the user. Init MUST define the following:
            1.) m_lambda    - the relaxation parameter
            2.) m_interpBC  - an array of InterpStencil objects. This is done using the helper function BaseOp::buildInterpStencils

            In addition, the user MAY define any other relevent stencils or temporaries which are needed for the implementation of the flux
            The user does NOT need to call init. It is called normally by define. 
        */
        virtual inline void init() = 0;

        /// Interpolation Stencil Helper
        /**
            A wrapper for Proto::InterpStencil<...>::Build. Creates the desired stencils for boundary value interpolation.
            The inputs to this function will likely need some tweaking after some user feedback. 

            \param a_maxShift       The maximum Manhatten distance for a cell to be included in the interpolation stencil's footprint (inclusive)
            \param a_footprint      The bounding Box for the interpolation stencil's footprint. 
            \param a_order          The desired polynomial order used for interpolation (e.g. 2 for quadratic). All terms are used by default.

            See the documentation for Proto::InterpStencil::Build for more information
        */
        inline void buildInterpStencils(int a_maxShift, Proto::Box a_footprint, int a_order);
        
        /// Interpolation Stencil Helper (Vector Version)
        /**
            Overload of buildInterpStencils which takes a std::vector of Proto::Point objects for the footprint

            \param a_footprint      std::vector of Proto::Point objects which defines the interpolation stencil's footprint.
            \param a_order          The desired polynomial order used for interpolation (e.g. 2 for quadratic). All terms are used by default.

            See the documentation for Proto::InterpStencil::Build for more information
        */
        inline void buildInterpStencils(std::vector<Proto::Point> a_footprint, int a_order);

        /// Get Number of Components
        static constexpr int numcomps() {return NCOMP;} 
        
        /// Get Number of Ghost Cells
        /**
            Refers to the number of ghost cells needed for the operator itself
        */
        static Proto::Point ghost() {return Proto::Point::Ones(NGHOST);}

        /// Get Number of Interpolation Ghost Cells
        /**
            Refers to the number of ghost cells needed for the boundary interpolation stencil.
        */
        inline Proto::Point interpGhost() {return this->m_interpBC[0].spanPoint();}
        
        /// Get Grid Spacing
        inline Real dx() const {return m_dx;}
        
        /// User Defined Flux Function
        /**
            Defines how the flux of the operator is calculated. This function - along with init - is what defines the core of any class derived from BaseOp. 
            The framework assumes flux is calculated ON THE LOWER FACE. That is to say, \f$F_{x,i}\equiv\frac{\partial \phi}{\partial x}_{x_{i-1/2}} \f$
            To facilitate this assumption, BaseOp provides the m_flux[DIM] stencils which implement the simple flux \f$F_{x,i}\approx\frac{\phi(x_i) - \phi(x_{i-1})}{dx}\f$.

            \param a_flux   A patch for the flux output.
            \param a_phi    A patch reference containing the independent variable data
            \param a_dir    Desired direction in which to evaluate the flux
            \param a_dx     Grid spacing to use when evaluating the flux
            \param a_scale  (Optional) value by which to scale the evaluated flux
        */
        virtual inline void flux(patch& a_flux, const patch& a_phi,
                int a_dir, Real a_dx, Real a_scale = 1.0) const = 0;
        /* OLD VERSION
        virtual inline patch flux(const patch& a_phi,
                int a_dir, Real a_dx, Real a_scale = 1.0) const = 0;
        */

        /// Apply Operator on Patch
        /**
            Apply the operator on a patch. The operator is evaluated as the divergence of the outputs generated by BaseOp::flux.
            
            \param a_rhs    A patch reference where the result is stored
            \param a_phi    A patch reference containing the independent variable data
            \param a_dx     Grid spacing used to evaluate the operator
            \param a_scale  (Optional) value by which to scale the evaluated operator
        */
        virtual inline void apply(patch& a_rhs, patch& a_phi,
                Real a_dx, Real a_scale = 1.0) const;
        
        //TODO: Make this private (user can edit the patch version)
        /// Apply Operator on Level
        /**
            Convenience function which calls the patch version of apply on each patch in a level
            
            Not for public use; modify the patch version in a derived class.

            \param a_rhs    A Chombo LevelData<FArrayBox> reference for the output    
            \param a_phi    A Chombo LevelData<FArrayBox> reference containing the independent variable data
            \param a_dx     Grid spacing used to evaluate the operator
            \param a_scale  (Optional) value by which to scale the evaluated operator
        */
        virtual inline void apply(LevelData<FArrayBox>& a_rhs,
                LevelData<FArrayBox>& a_phi,
                Real a_dx, Real a_scale = 1.0) const;
        
        //TODO: Should be private
        /// Add Operator To Level
        /**
            Convenience function which applies the operator and adds the result to an existing patch for each patch in a level.
            
            Not for public use; modify the patch version in a derived class.
            
            \param a_rhs    A Chombo LevelData<FArrayBox> reference to which the evaluated operator will be added
            \param a_phi    A Chombo LevelData<FArrayBox> reference containing the independent variable data
            \param a_dx     Grid spacing used to evaluate the operator
            \param a_scale  (Optional) value by which to scale the evaluated operator (does not scale the input data, only the increment)
        */
        virtual inline void plusApply(LevelData<FArrayBox>& a_rhs,
                LevelData<FArrayBox>& a_phi,
                Real a_dx, Real a_scale = 1.0) const;
        
        /// Compute Residual on Patch
        /**
            Computes the residual on a patch as \f$ r \equiv R - L(\phi) \f$ where r is the residual, R is the right-hand size, and L is the operator.
            Stores the result in a given patch.
            
            \param a_res    A patch reference containing the output residual data
            \param a_phi    A patch reference containing the independent variable data
            \param a_R      A patch reference containing the right-hand size data
            \param a_dx     Grid spacing used to evaluate the operator
        */
        virtual inline void residual(patch& a_res,
                patch& a_phi, const patch& a_R, const Real& a_dx) const;
        
        //TODO: Should be private
        /// Compute Residual Max Norm on Level
        /**
            Convenience function that computes the residual on each patch in a level and stores the output in a given LevelData.
            Outputs the max norm of all residuals on the level. 

            Not for public use; modify the patch version in a derived class.
            
            \param a_res    A Chombo LevelData<FArrayBox> reference containing the output residual data
            \param a_phi    A Chombo LevelData<FArrayBox> reference containing the independent variable data
            \param a_R      A Chombo LevelData<FArrayBox> reference containing the right-hand size data
            \param a_dx     Grid spacing used to evaluate the operator
        */
        virtual inline Real residual(LevelData<FArrayBox>& a_res,
                LevelData<FArrayBox>& a_phi,
                const LevelData<FArrayBox>& a_R,
                const Real a_dx) const;
        
        //TODO: Should be private
        /// Compute Residual Max Norm on Level
        /**
            Just computes the residual max-norm on a level without storing the residual field

            Not for public use; modify the patch version in a derived class.
            
            \param a_phi    A Chombo LevelData<FArrayBox> reference containing the independent variable data
            \param a_R      A Chombo LevelData<FArrayBox> reference containing the right-hand size data
            \param a_dx     Grid spacing used to evaluate the operator
        */
        virtual inline Real residual(LevelData<FArrayBox>& a_phi,
                const LevelData<FArrayBox>& a_R,
                Real a_dx) const;

        //TODO: Should be private, possibly removed
        /// Compute Residual on Level
        /**
            Convenience overload which uses the default grid spacing

            Not for public use; modify the patch version in a derived class.
            
            \param a_res    A Chombo LevelData<FArrayBox> reference containing the output residual data
            \param a_phi    A Chombo LevelData<FArrayBox> reference containing the independent variable data
            \param a_R      A Chombo LevelData<FArrayBox> reference containing the right-hand size data
        */
        virtual inline Real residual(LevelData<FArrayBox>& a_res,
                LevelData<FArrayBox>& a_phi,
                const LevelData<FArrayBox>& a_R) const;
        
        //TODO: Should be private, possibly removed
        /// Compute Residual Max Norm on Level
        /**
            Convenience overload which uses the default grid spacing
            
            Not for public use; modify the patch version in a derived class.

            \param a_phi    A Chombo LevelData<FArrayBox> reference containing the independent variable data
            \param a_R      A Chombo LevelData<FArrayBox> reference containing the right-hand size data
        */
        virtual inline Real residual(LevelData<FArrayBox>& a_phi,
                const LevelData<FArrayBox>& a_R) const;

        //TODO: Should be private with a public patch version?
        /// Level Relax
        /**
            Relax independent variable data using Jacobi Iteration.
            The user may override this function to implement their own iterative scheme (e.g. GSRB)
            
            \param a_phi    Chombo LevelData<FArrayBox> reference containing the independent variable data to be relaxed
            \param a_R      Chombo LevelData<FArrayBox> reference containing the right-hand size data used for relaxation
            \param a_n      Number of iterations
        */
        virtual inline void relax(LevelData<FArrayBox>& a_phi,
                const LevelData<FArrayBox>& a_R,
                int a_n) const;
        
        //TODO: Should be private with a public patch version?
        /// Level Coarsen / Restrict
        /**
            Coarsens input data using a conservative average: \f$ \phi_C \equiv \frac{1}{N}\Sigma_{i = 1}^{N}(\phi) \f$ where \f$ N \equiv refRatio^{DIM} \f$
            The fine data is assumed to be on a THIS level, and to be coarsened onto a level coarser by either AMR_REFRATIO (if *this was defined with amr = true) or MG_REFRATIO

            \param a_phiC   Chombo LevelData<FArrayBox> reference to the coarsened data
            \param a_phi    Chombo LevelData<FArrayBox> reference to the fine data
        */
        virtual inline void coarsen(LevelData<FArrayBox>& a_phiC,
                const LevelData<FArrayBox>& a_phi);
        
        //TODO: Should be private
        /// Multigrid Coarse Right-hand Side
        /**
            Convenience function that computes the right-hand side for the next coarser level for vanilla multigrid.
            \f$ R_C \equiv L(\phi_C) + \langle R - L(\phi)\rangle\f$.

            \param  a_RC    Chombo LevelData<FArrayBox> reference to the computed coarse right-hand side
            \param  a_phiC  Chombo LevelData<FArrayBox> reference to the coarse independent variable data
            \param  a_phi   Chombo LevelData<FArrayBox> reference to the independent variable data on this level
            \param  a_R     Chombo LevelData<FArrayBox> reference to the right-hand side data on this level
        */
        virtual inline void coarseRhs(LevelData<FArrayBox>& a_RC,
                LevelData<FArrayBox>& a_phiC,
                LevelData<FArrayBox>& a_phi,
                const LevelData<FArrayBox>& a_R);
        
        //TODO: Should be private
        /// Fine Correction
        /**
            Updates the independent variable data on this level using stored and relaxed coarse data: \f$ \phi += \mathcal{I}(\phi_C - \phi_{C,0}) \f$

            \param  a_phi   Chombo LevelData<FArrayBox> reference to the independent variable data on this level
            \param  a_phiC  Chombo LevelData<FArrayBox> reference to the relaxed coarse independent variable data
            \param  a_phiC0 Chombo LevelData<FArrayBox> reference to the stored independent variable data
        */
        virtual inline void fineCorrection(LevelData<FArrayBox>& a_phi,
                LevelData<FArrayBox>& a_phiC,
                const LevelData<FArrayBox>& a_phiC0);

        //TODO: Should be private
        /// Interpolate Boundary
        /**
            Interpolate values from coarse data to fine data ghost cells. The order of accuracy of interpolation is based on the inputs to buildInterpStencils which
            the user should call in their implementation of init(). If the user does not initialize these stencils themselves, the interpolation will default to
            3rd order method.

            \param  a_phi   Chombo LevelData<FArrayBox> reference to the independent variable data on this level to which data will be interpolated
            \param  a_phiC  Chombo LevelData<FArrayBox> reference to the independent variable data on the next coarser level from which data will be interplated
        */
        virtual inline void interpBoundary(LevelData<FArrayBox>& a_phi,
                const LevelData<FArrayBox>& a_phiC);
        
        //TODO: Should be private
        /// Interpolate Multigrid Boundary
        /**
            Version of interpolateBoundary used to compute boundary information on Multigrid levels between AMR levels. Not recommended for public use.
            
            \param  a_phi   Chombo LevelData<FArrayBox> reference to the independent variable data on this level to which data will be interpolated
            \param  a_phiC  Chombo LevelData<FArrayBox> reference to the independent variable data on the next coarser level from which data will be interplated
            \param  a_level Multigrid level on to which data is being interpolated. A value between 0 and log2(AMR_REFRATIO/MG_REFRATIO);
        */
        virtual inline void interpBoundary(LevelData<FArrayBox>& a_phi,
                LevelData<FArrayBox>& a_temp, int a_level);
        
        //TODO: Should be private
        /// Boundary Reflux
        /**
            This function is used internally by AMRFAS and it is not recommended for public use. 

            Implements the correction: \f$ R_C = R_C + \nabla \cdot \delta F \mid \delta\mathcal{C}(\Omega_f) \f$ where \f$ \delta F \equiv \langle F(\phi) \rangle - F(\phi_C) \f$ 
            This correction has the effect of fixing the flux of any terms of the form \f$ L(\phi_C) \f$ in \f$ R_C \f$ so as to be consistent with the flux at the next finer level (e.g. this level).
            
            Note that the data in \f$ \mathcal{C}(\Omega_f)\f$  will be corrupted by this function and should be overwritten or not used. 

            \param  a_RC    Chombo LevelData<FArrayBox> reference to the coarse data being corrected
            \param  a_phiC  Chombo LevelData<FArrayBox> reference to the independent variable data on the next coarser level
            \param  a_phi   Chombo LevelData<FArrayBox> reference to the independent variable data on this level
            \param  a_flux  A Chombo LevelFluxRegister object built from the DisjointBoxLayouts on this and the next coarser level
            \param  a_scale (Optional) Scales the correction. Mostly used to flip the sign (e.g. when correcting a residual quantity)

        */
        virtual inline void reflux(LevelData<FArrayBox>& a_RC,
                LevelData<FArrayBox>& a_phiC,
                LevelData<FArrayBox>& a_phi,
                LevelFluxRegister& a_flux,
                Real a_scale = 1.0);
        
        //TODO: Should be private
        /// AMR Coarse Residual
        /**
            Computes the residual on the next coarser level of an AMR hierarchy including boundary refluxing. 

            \param a_resC   Chombo LevelData<FArrayBox> reference to the coarse residual
            \param a_rhoC   Chombo LevelData<FArrayBox> reference to the coarse right-hand side
            \param a_phiC   Chombo LevelData<FArrayBox> reference to the coarse independent variable data
            \param a_resC   Chombo LevelData<FArrayBox> reference to the right-hand side on this level
            \param a_resC   Chombo LevelData<FArrayBox> reference to the independent variable data on this level
            \param  a_flux  A Chombo LevelFluxRegister object built from the DisjointBoxLayouts on this and the next coarser level
        */
        virtual inline void coarseResidual(LevelData<FArrayBox>& a_resC,
                const LevelData<FArrayBox>& a_rhoC,
                LevelData<FArrayBox>& a_phiC,
                const LevelData<FArrayBox>& a_R,
                LevelData<FArrayBox>& a_phi,
                LevelFluxRegister& a_flux);
    protected:
        
        Proto::Stencil<Real>        m_flux[DIM];    ///< Simple stencil computing the undivided difference \f$ \phi_{i-1/2}\equiv \phi_i - \phi_{i-1} \f$
        Proto::Stencil<Real>        m_div[DIM];     ///< Simple stencil computing the undivided difference \f$ \phi_{i+1/2}\equiv \phi_{i+1} - \phi_i \f$
        Proto::Stencil<Real>        m_avg;          ///< Stencil for computing a conservative average
        Proto::InterpStencil<Real>  m_interp;       ///< Piecewise constant InterStencil. 
        std::vector<Proto::InterpStencil<Real>>  m_interpBC;    ///< Set of InterpStencil used for interpolating boundary conditions at different ratios
        
        Real                        m_dx;           ///< Grid spacing for this level
        Real                        m_lambda;       ///< relaxation parameter
        unsigned int                m_refRatio;     ///< refinement ratio between this level and the next coarser level
        LevelData<FArrayBox>             m_coarseTemp;   ///< temporary data structure used for regridding between this level and the next coarser level
        bool                        m_defined;      ///< true if define(...) has been called
};

template<unsigned int NCOMP, unsigned int NGHOST>
void BaseOp<NCOMP, NGHOST>::define(const DisjointBoxLayout& a_layout,
    Real a_dx, bool a_amr)
{
    PR_TIME("BaseOp::define");
    m_dx = a_dx;
    if (a_amr)
    {
        m_refRatio = AMR_REFRATIO;
    } else {
        m_refRatio = MG_REFRATIO;
    }
    for (int ii = 0; ii < DIM; ii++)
    {
        m_flux[ii] = 1.0*Proto::Shift::Zeros()
                   - 1.0*Proto::Shift::Basis(ii,-1);
        m_div[ii] = 1.0*Proto::Shift::Basis(ii,1)
                  - 1.0*Proto::Shift::Zeros();
    }
    m_avg = Proto::Stencil<Real>::AvgDown(m_refRatio);
    m_interp = Proto::InterpStencil<Real>::PiecewiseConstant(Proto::Point::Ones(m_refRatio));
    m_lambda = 0;
    init(); //define lambda and any other custom quantities.
    if (m_interpBC.size() == 0)
    {
        // default implementation if none is provided through init()
        this->buildInterpStencils(DIM, Proto::Box::Cube(3).shift(Proto::Point::Ones(-1)), 2);
    }
    assert(m_lambda != 0);
    
    Proto::Point tempGhost;
    if (a_amr)
    {
        tempGhost = this->ghost() + this->interpGhost();
    } else {
        tempGhost = this->ghost();
    }
    DisjointBoxLayout coarseLayout;
    coarsen_dbl(coarseLayout, a_layout, m_refRatio);
    m_coarseTemp.define(coarseLayout, NCOMP, tempGhost);
    m_defined = true;
}

template<unsigned int NCOMP, unsigned int NGHOST>
void BaseOp<NCOMP, NGHOST>::buildInterpStencils(int a_maxShift, Proto::Box a_footprint, int a_order)
{
    int r = MG_REFRATIO;
    while (r <= AMR_REFRATIO)
    {
        m_interpBC.push_back(Proto::InterpStencil<Real>::Build(a_maxShift, a_footprint, a_order, r));
        r *= MG_REFRATIO;
    }
}

template<unsigned int NCOMP, unsigned int NGHOST>
void BaseOp<NCOMP, NGHOST>::buildInterpStencils(std::vector<Proto::Point> a_footprint, int a_order)
{
    std::remove(a_footprint.begin(), a_footprint.end(), Proto::Point::Zeros());
    int r = MG_REFRATIO;
    while (r <= AMR_REFRATIO)
    {
        m_interpBC.push_back(Proto::InterpStencil<Real>::Build(a_footprint, a_order, r));
        r *= MG_REFRATIO;
    }
}

template<unsigned int NCOMP, unsigned int NGHOST>
void BaseOp<NCOMP, NGHOST>::apply(patch& a_rhs, patch& a_phi,
        Real a_dx, Real a_scale) const
{
    PR_TIME("BaseOp::apply (on a patch)");
    PROTO_ASSERT(m_defined, "Operator is not defined!");
    for (int ii = 0; ii < DIM; ii++)
    {
        patch flx; //should make this a member of BaseOp
        this->flux(flx, a_phi, ii, a_dx, 1.0);
        if (ii == 0)
        {
            a_rhs |= this->m_div[ii](flx, a_scale/a_dx);
        } else {
            a_rhs += this->m_div[ii](flx, a_scale/a_dx);
        }
    }
}

template<unsigned int NCOMP,  unsigned int NGHOST>
void BaseOp<NCOMP, NGHOST>::apply(LevelData<FArrayBox>& a_rhs,
        LevelData<FArrayBox>& a_phi,
        Real a_dx, Real a_scale) const
{
    PR_TIME("BaseOp::apply (on a level)");
    a_phi.exchange();
    auto iter = a_phi.dataIterator();
    for (iter.begin(); iter.ok(); ++iter)
    {
        //patch phi = a_phi[iter];
        //patch rhs = a_rhs[iter];
        patch phi, rhs;
        ProtoCh::aliasBoxData<Real, NCOMP>(phi, a_phi[iter]);
        ProtoCh::aliasBoxData<Real, NCOMP>(rhs, a_rhs[iter]);

        apply(rhs, phi, a_dx, a_scale);
    }
}

template<unsigned int NCOMP,  unsigned int NGHOST>
void BaseOp<NCOMP, NGHOST>::plusApply(LevelData<FArrayBox>& a_rhs,
        LevelData<FArrayBox>& a_phi,
        Real a_dx, Real a_scale) const
{
    PR_TIME("BaseOp::plusApply");
    a_phi.exchange();
    auto layout = a_phi.disjointBoxLayout();
    auto iter = a_phi.dataIterator();
    for (iter.begin(); iter.ok(); ++iter)
    {
        //patch phi = a_phi[iter];
        //patch rhs = a_rhs[iter];
        patch phi, rhs;
        ProtoCh::aliasBoxData<Real, NCOMP>(phi, a_phi[iter]);
        ProtoCh::aliasBoxData<Real, NCOMP>(rhs, a_rhs[iter]);

        patch temp(layout[iter]);
        apply(temp, phi, a_dx, a_scale);
        rhs += temp;
    }
}

template<unsigned int NCOMP, unsigned int NGHOST>
void BaseOp<NCOMP, NGHOST>::residual(patch& a_res,
        patch& a_phi,
        const patch& a_R,
        const Real& a_dx) const
{
    apply(a_res, a_phi, a_dx, -1.0);
    a_res += a_R;
}

template<unsigned int NCOMP, unsigned int NGHOST>
Real BaseOp<NCOMP, NGHOST>::residual(LevelData<FArrayBox>& a_res,
        LevelData<FArrayBox>& a_phi,
        const LevelData<FArrayBox>& a_R,
        Real a_dx) const
{
    PR_TIME("BaseOp::residual");
    a_phi.exchange();
    Real resnorm = 0.0;
    auto iter = a_phi.dataIterator();
    for (iter.begin(); iter.ok(); ++iter)
    {
        //patch phi = a_phi[iter];
        //patch rhs = a_R[iter];
        //patch res = a_res[iter];
        patch phi, rhs, res;

        ProtoCh::aliasBoxData<Real, NCOMP>(phi, a_phi[iter]);
        ProtoCh::aliasBoxData<Real, NCOMP>(rhs, a_R[iter]);
        ProtoCh::aliasBoxData<Real, NCOMP>(res, a_res[iter]);
        
        residual(res, phi, rhs, a_dx);
        resnorm = std::max(res.absMax(), resnorm);
    }
    return resnorm;
}

template<unsigned int NCOMP, unsigned int NGHOST>
Real BaseOp<NCOMP, NGHOST>::residual(LevelData<FArrayBox>& a_phi,
        const LevelData<FArrayBox>& a_R,
        Real a_dx) const
{
    PR_TIME("BaseOp::residual");
    a_phi.exchange();
    Real resnorm = 0.0;
    auto iter = a_phi.dataIterator();
    for (iter.begin(); iter.ok(); ++iter)
    {
        //patch phi = a_phi[iter];
        //patch R = a_R[iter];
        patch phi, R;
        ProtoCh::aliasBoxData<Real, NCOMP>(phi, a_phi[iter]);
        ProtoCh::aliasBoxData<Real, NCOMP>(R, a_R[iter]);

        patch res(R.box());
        residual(res, phi, R, a_dx);
        resnorm = std::max(res.absMax(), resnorm);
    }
    return resnorm;
}

template<unsigned int NCOMP, unsigned int NGHOST>
Real BaseOp<NCOMP, NGHOST>::residual(LevelData<FArrayBox>& a_res,
        LevelData<FArrayBox>& a_phi,
        const LevelData<FArrayBox>& a_R ) const
{
    return residual(a_res, a_phi, a_R, m_dx);
}

template<unsigned int NCOMP, unsigned int NGHOST>
Real BaseOp<NCOMP, NGHOST>::residual(LevelData<FArrayBox>& a_phi,
        const LevelData<FArrayBox>& a_R ) const
{
    return residual(a_phi, a_R, m_dx);
}

template<unsigned int NCOMP, unsigned int NGHOST>
void BaseOp<NCOMP, NGHOST>::relax(LevelData<FArrayBox>& a_phi,
        const LevelData<FArrayBox>& a_R,
        int a_n) const
{
    PR_TIME("BaseOp::relax");
    
    auto layout = a_phi.disjointBoxLayout();
    auto iter = a_phi.dataIterator();
    for (int ii = 0; ii < a_n; ii++)
    {
        //writeLevel("Phi_relax_%i.hdf5", ii);
        a_phi.exchange();
        for (iter.begin(); iter.ok(); ++iter)
        {
            //patch phi = a_phi[iter];
            //patch R = a_R[iter];
            patch phi, R;
            ProtoCh::aliasBoxData<Real, NCOMP>(phi, a_phi[iter]);
            ProtoCh::aliasBoxData<Real, NCOMP>(R, a_R[iter]);
            patch res(layout[iter]);
           
            residual(res, phi, R, m_dx); 
            res *= m_lambda;
            phi += res;
        }
    }
}

template<unsigned int NCOMP, unsigned int NGHOST>
void BaseOp<NCOMP, NGHOST>::coarsen(LevelData<FArrayBox>& a_phiC,
        const LevelData<FArrayBox>& a_phi)
{
    PR_TIME("BaseOp::coarsen");
    PROTO_ASSERT(m_defined, "Operator is not defined!");
    auto iter = a_phi.dataIterator();
    for (iter.begin(); iter.ok(); ++iter)
    {
        //patch phi = a_phi[iter];
        //patch temp(m_coarseTemp[iter]);
        patch phi, temp;
        ProtoCh::aliasBoxData<Real, NCOMP>(phi, a_phi[iter]);
        ProtoCh::aliasBoxData<Real, NCOMP>(temp, m_coarseTemp[iter]);
        
        temp |= m_avg(phi);
    }
    m_coarseTemp.copyTo(a_phiC);
}

template<unsigned int NCOMP, unsigned int NGHOST>
void BaseOp<NCOMP, NGHOST>::coarseRhs(LevelData<FArrayBox>& a_RC,
        LevelData<FArrayBox>& a_phiC,
        LevelData<FArrayBox>& a_phi,
        const LevelData<FArrayBox>& a_R)
{
    PR_TIME("BaseOp::coarseRhs");
    a_phi.exchange();
    a_phiC.exchange();
    
    auto fineLayout = a_phi.disjointBoxLayout();
    auto fiter = a_phi.dataIterator();
    for (fiter.begin(); fiter.ok(); ++fiter)
    {
        //patch phi = a_phi[fiter];
        //patch R = a_R[fiter];
        //patch temp(m_coarseTemp[fiter]);
        patch phi, R, temp;
        ProtoCh::aliasBoxData<Real, NCOMP>(phi, a_phi[fiter]);
        ProtoCh::aliasBoxData<Real, NCOMP>(R, a_R[fiter]);
        ProtoCh::aliasBoxData<Real, NCOMP>(temp, m_coarseTemp[fiter]);

        patch res(fineLayout[fiter]);
        
        residual(res, phi, R, m_dx);
        temp |= m_avg(res);
    }
    m_coarseTemp.copyTo(a_RC);
    plusApply(a_RC, a_phiC, m_dx*m_refRatio);
}

template<unsigned int NCOMP, unsigned int NGHOST>
void BaseOp<NCOMP, NGHOST>::fineCorrection(LevelData<FArrayBox>& a_phi,
        LevelData<FArrayBox>& a_phiC,
        const LevelData<FArrayBox>& a_phiC0)
{
    PR_TIME("BaseOp::fineCorrection");
    PROTO_ASSERT(m_defined, "Operator is not defined!");
    auto citer = a_phiC.dataIterator();
    for (citer.begin(); citer.ok(); ++citer)
    {
        //patch phiC = a_phiC[citer];
        //patch phiC0 = a_phiC0[citer];
        patch phiC, phiC0;
        ProtoCh::aliasBoxData<Real, NCOMP>(phiC, a_phiC[citer]);
        ProtoCh::aliasBoxData<Real, NCOMP>(phiC0, a_phiC0[citer]);
        
        phiC0 -= phiC;
        phiC0 *= (-1.0);
    }
    a_phiC0.copyTo(m_coarseTemp);
    auto fineLayout = a_phi.disjointBoxLayout();
    auto fiter = a_phi.dataIterator();
    for (fiter.begin(); fiter.ok(); ++fiter)
    {
        //patch phi = a_phi[fiter];
        //patch temp(m_coarseTemp[fiter]);
        patch phi, temp;
        ProtoCh::aliasBoxData<Real, NCOMP>(phi, a_phi[fiter]);
        ProtoCh::aliasBoxData<Real, NCOMP>(temp, m_coarseTemp[fiter]);
        
        Proto::Box b = Proto::Box(fineLayout[fiter]).coarsen(m_refRatio);
        phi += m_interp(temp, b);
    }
}

template<unsigned int NCOMP, unsigned int NGHOST>
void BaseOp<NCOMP, NGHOST>::interpBoundary(LevelData<FArrayBox>& a_phi,
        const LevelData<FArrayBox>& a_phiC)
{
    a_phiC.copyTo(m_coarseTemp);
    interpBoundary(a_phi, m_coarseTemp, m_interpBC.size() - 1);
}

template<unsigned int NCOMP, unsigned int NGHOST>
void BaseOp<NCOMP, NGHOST>::interpBoundary(LevelData<FArrayBox>& a_phi,
        LevelData<FArrayBox>& a_temp, int a_level)
{
    PR_TIME("BaseOp::interpBoundary");
    PROTO_ASSERT(m_defined, "Operator is not defined!");
    Proto::InterpStencil<Real>& interpStencil = m_interpBC[a_level];
    auto iter = a_phi.dataIterator();
    auto fineLayout = a_phi.disjointBoxLayout();
    Proto::Box fineDomain = fineLayout.physDomain().domainBox();

    int r[DIM];
    for (int ii = 0; ii < DIM; ii++)
    {
        r[ii] = std::min((int)fineDomain.size(ii), MAXBOXSIZE);
    }
    Proto::Point bitRatio(r);
    Proto::Box bitBox = fineDomain.coarsen(bitRatio);

    Proto::Box U = Proto::Box::Cube(3).shift(Proto::Point::Ones(-1));
    for (iter.begin(); iter.ok(); ++iter)
    {
        //patch coarseData = a_temp[iter];
        //patch fineData = a_phi[iter];
        patch coarseData, fineData;
        ProtoCh::aliasBoxData<Real, NCOMP>(coarseData, a_temp[iter]); 
        ProtoCh::aliasBoxData<Real, NCOMP>(fineData, a_phi[iter]); 
        
        Proto::Box coarsePatch = coarseData.box();
        Proto::Point bitPoint = Proto::Box(fineLayout[iter]).low()/bitRatio;
        for (auto niter = U.begin(); niter != U.end(); ++niter)
        {
            if (*niter == Proto::Point::Zeros()) {continue;}
            Proto::Point bitNeighbor = bitPoint + *niter;
            if (!bitBox.contains(bitNeighbor))
            {
                Proto::Box coarseNeighborBox(bitNeighbor, bitNeighbor);
                coarseNeighborBox =
                        coarseNeighborBox.refine(bitRatio/m_refRatio);
                Proto::Box intersect = coarseNeighborBox & coarsePatch;
                fineData |= interpStencil(coarseData, intersect);
            }
        }
    }
}

template<unsigned int NCOMP, unsigned int NGHOST>
void BaseOp<NCOMP, NGHOST>::reflux(LevelData<FArrayBox>& a_RC,
        LevelData<FArrayBox>& a_phiC,
        LevelData<FArrayBox>& a_phi,
        LevelFluxRegister& a_flux,
        Real a_scale)
{
    PR_TIME("BaseOp::reflux");
    //Turn this on to produce debugging outputs from reflux
    //bool debug = false;
    PROTO_ASSERT(m_defined, "Operator is not defined!");
    
    a_flux.setToZero();
    a_phi.exchange();
    a_phiC.exchange();

    Interval comps(0, NCOMP-1);
   
    auto fineLayout     = a_phi.disjointBoxLayout();
    auto coarseLayout   = a_phiC.disjointBoxLayout();
    auto fiter = a_phi.dataIterator();
    auto citer = a_phiC.dataIterator();
    Real cdx = m_dx*m_refRatio;
    Proto::Box invalidRegion = fineLayout.physDomain().domainBox();
    invalidRegion = invalidRegion.coarsen(m_refRatio);
    
    for (fiter.begin(); fiter.ok(); ++fiter)
    {
        Proto::Box patchBox = fineLayout[fiter];
        //patch phi = a_phi[fiter];
        patch phi;
        ProtoCh::aliasBoxData<Real, NCOMP>(phi, a_phi[fiter]);
        for (int dir = 0; dir < DIM; dir++)
        {
            patch fineFlux;
            this->flux(fineFlux, phi, dir, m_dx);
            FArrayBox fineFluxFAB(fineFlux, Proto::Point::Basis(dir));
            a_flux.incrementFine(fineFluxFAB, 1.0, fiter(), comps, comps, dir);
        }
    }
    for (citer.begin(); citer.ok(); ++citer)
    {
        Proto::Box patchBox = coarseLayout[citer];
        //patch phiC = a_phiC[citer];
        patch phiC;
        ProtoCh::aliasBoxData<Real, NCOMP>(phiC, a_phiC[citer]);
        for (int dir = 0; dir < DIM; dir++)
        {
            patch coarseFlux;
            this->flux(coarseFlux,phiC, dir, cdx);
            FArrayBox coarseFluxFAB(coarseFlux, Proto::Point::Basis(dir));
            a_flux.incrementCoarse(coarseFluxFAB, 1.0, citer(), comps, comps, dir);
        }
    }
    // negative sign -> defaults to refluxing with <F_fine> - F_coarse
    // yes, I am quite sure that this is the way LFR is coded -CLG
    Real scale = -a_scale/cdx;
    a_flux.reflux(a_RC, scale);
    /*
    if (debug)
    {
        LevelData<FArrayBox> LPhi, LPhiC, Flux, FluxC_X, FluxC_Y, FluxF_X, FluxF_Y;
        LPhi.define(fineLayout, numcomps(), Proto::Point::Zeros());
        LPhiC.define(coarseLayout, numcomps(), Proto::Point::Zeros());
        Flux.define(coarseLayout, numcomps(), Proto::Point::Zeros());
        FluxC_X.define(coarseLayout, numcomps(), Proto::Point::Ones());
        FluxC_Y.define(coarseLayout, numcomps(), Proto::Point::Ones());
        FluxF_X.define(fineLayout, numcomps(), Proto::Point::Ones());
        FluxF_Y.define(fineLayout, numcomps(), Proto::Point::Ones());
        for (citer.reset(); citer.ok(); ++citer)
        {
            patch lphi = LPhiC[citer];
            patch phi = a_phiC[citer];
            patch flx = Flux[citer];
            patch flx_x = FluxC_X[citer];
            patch flx_y = FluxC_Y[citer];
            flx.setVal(0);
            flx_x.setVal(0);
            flx_y.setVal(0);
            auto tmp_x = this->flux(phi, 0, cdx);
            tmp_x.copyTo(flx_x);
            auto tmp_y = this->flux(phi, 1, cdx);
            tmp_y.copyTo(flx_y);
            apply(lphi, phi, cdx, 1.0);
        }
        for (fiter.reset(); fiter.ok(); ++fiter)
        {
            patch lphi = LPhi[fiter];
            patch phi = a_phi[fiter];
            patch flx_x = FluxF_X[fiter];
            patch flx_y = FluxF_Y[fiter];
            flx_x.setVal(0);
            flx_y.setVal(0);
            auto tmp_x = this->flux(phi, 0, m_dx);
            tmp_x.copyTo(flx_x);
            auto tmp_y = this->flux(phi, 1, m_dx);
            tmp_y.copyTo(flx_y);

            apply(lphi, phi, m_dx, 1.0);
        }
        LPhiC.copyTo(m_coarseTemp);
        a_flux.reflux(Flux, scale);
        Real intFine = integrate(LPhi, m_dx);
        Real intCoarse = integrate(m_coarseTemp, cdx);
        std::cout << "Coarse integral: " << intCoarse << std::endl;
        std::cout << "Fine integral: " << intFine << std::endl;
        std::cout << "Fine - Coarse: " << intFine - intCoarse << std::endl;
        std::cout << "Flux register integral: " << integrate(Flux, cdx) << std::endl;
        writeLevel(Flux, "reflux_register.hdf5");
        writeLevel(FluxC_X, "reflux_coarseFlux_X.hdf5");
        writeLevel(FluxC_Y, "reflux_coarseFlux_Y.hdf5");
        writeLevel(FluxF_X, "reflux_fineFlux_X.hdf5");
        writeLevel(FluxF_Y, "reflux_fineFlux_Y.hdf5");
    }
    */
}

template<unsigned int NCOMP, unsigned int NGHOST>
void BaseOp<NCOMP, NGHOST>::coarseResidual(LevelData<FArrayBox>& a_resC,
        const LevelData<FArrayBox>& a_rhoC,
        LevelData<FArrayBox>& a_phiC,
        const LevelData<FArrayBox>& a_R,
        LevelData<FArrayBox>& a_phi,
        LevelFluxRegister& a_flux)
{
    PR_TIME("BaseOp::coarseResidual");
    Real cdx = m_dx*m_refRatio;
    coarsen(a_phiC, a_phi);
    interpBoundary(a_phi, a_phiC);
    a_phiC.exchange();
    a_phi.exchange();
    
    auto coarseLayout = a_resC.disjointBoxLayout();
    auto citer = a_resC.dataIterator();
    for (citer.begin(); citer.ok(); ++citer)
    {
        //patch resC = a_resC[citer];
        //patch rhoC = a_rhoC[citer];
        //patch phiC = a_phiC[citer];
        patch resC, rhoC, phiC;
        ProtoCh::aliasBoxData<Real, NCOMP>(resC, a_resC[citer]);
        ProtoCh::aliasBoxData<Real, NCOMP>(rhoC, a_rhoC[citer]);
        ProtoCh::aliasBoxData<Real, NCOMP>(phiC, a_phiC[citer]);

        apply(resC, phiC, cdx, -1.0);
        resC += rhoC; 
    }
    reflux(a_resC, a_phiC, a_phi, a_flux,1.0);
    
    auto fineLayout = a_phi.disjointBoxLayout();
    auto fiter = a_phi.dataIterator();
    for (fiter.begin(); fiter.ok(); ++fiter)
    {
        //patch temp = m_coarseTemp[fiter];
        //patch phi = a_phi[fiter];
        //patch R = a_R[fiter];
        patch temp, phi, R;
        ProtoCh::aliasBoxData<Real, NCOMP>(phi, a_phi[fiter]);
        ProtoCh::aliasBoxData<Real, NCOMP>(R, a_R[fiter]);
        ProtoCh::aliasBoxData<Real, NCOMP>(temp, m_coarseTemp[fiter]);
        
        patch res(fineLayout[fiter]);
        temp.setVal(0);

        residual(res, phi, R, m_dx);
        temp |= m_avg(res);
    }
    m_coarseTemp.copyTo(a_resC);
}
