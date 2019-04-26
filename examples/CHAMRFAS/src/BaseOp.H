#pragma once

#include "Proto.H"
#include "LevelData.H"
#include "AMRUtils.H"
#include "DisjointBoxLayout.H"
#include "LevelFluxRegister.H"
#include "Proto_WriteBoxData.H"
#include <cmath>

template<unsigned int NCOMP, unsigned int NGHOST, class DATA>
class BaseOp
{
    public:

        typedef Proto::BoxData<Real, NCOMP> patch;
        typedef Proto::Var<Real, NCOMP> var;

        // no constructors are defined; use weak construction. 
        inline BaseOp(){}
        inline void define(const DisjointBoxLayout& a_layout,
                Real a_dx, bool a_amr = true);
        
        virtual inline void init() = 0;

        static constexpr int numcomps() {return NCOMP;} 
        static constexpr Proto::Point ghost() {return Proto::Point::Ones(NGHOST);}
        inline Real dx() const {return m_dx;}
        virtual inline patch flux(const patch& a_phi,
                int a_dir, Real a_dx, Real a_scale = 1.0) const = 0;
        virtual inline void apply(patch& a_rhs, patch& a_phi,
                Real a_dx, Real a_scale = 1.0) const;
        virtual inline void apply(LevelData<DATA>& a_rhs,
                LevelData<DATA>& a_phi,
                Real a_dx, Real a_scale = 1.0) const;

        virtual inline Real residual(LevelData<DATA>& a_res,
            LevelData<DATA>& a_phi,
            const LevelData<DATA>& a_R,
            Real a_dx) const;
        virtual inline Real residual(LevelData<DATA>& a_phi,
            const LevelData<DATA>& a_R,
            Real a_dx) const;
        virtual inline Real residual(LevelData<DATA>& a_res,
            LevelData<DATA>& a_phi,
            const LevelData<DATA>& a_R) const;
        virtual inline Real residual(LevelData<DATA>& a_phi,
            const LevelData<DATA>& a_R) const;
       
        virtual inline void relax(LevelData<DATA>& a_phi,
            const LevelData<DATA>& a_R,
            int a_n) const;
        
        virtual inline void coarsen(LevelData<DATA>& a_phiC,
            const LevelData<DATA>& a_phi);
        virtual inline void coarseRhs(LevelData<DATA>& a_RC,
            LevelData<DATA>& a_phiC,
            LevelData<DATA>& a_phi,
            const LevelData<DATA>& a_R);
        virtual inline void fineCorrection(LevelData<DATA>& a_phi,
            LevelData<DATA>& a_phiC,
            const LevelData<DATA>& a_phiC0); 
            
        // AMR functions
        virtual inline void interpBoundary(LevelData<DATA>& a_phi,
            const LevelData<DATA>& a_phiC);
        virtual inline void interpBoundary(LevelData<DATA>& a_phi,
            LevelData<DATA>& a_temp,
            const Proto::InterpStencil<Real>& a_interp);
        virtual inline void coarseRhs(LevelData<DATA>& a_RC,
            LevelData<DATA>& a_phiC,
            LevelData<DATA>& a_phi,
            const LevelData<DATA>& a_R,
            const LevelData<DATA>& a_rhoC,
            LevelFluxRegister& a_flux);
        virtual inline void reflux(LevelData<DATA>& a_RC,
            LevelData<DATA>& a_phiC,
            LevelData<DATA>& a_phi,
            LevelFluxRegister& a_flux,
            Real a_scale = -1.0);
        virtual inline void refluxUndivided(LevelData<DATA>& a_RC,
            LevelData<DATA>& a_phiC,
            LevelData<DATA>& a_phi,
            LevelFluxRegister& a_flux,
            Real a_scale = -1.0);
        virtual inline void coarseResidual(LevelData<DATA>& a_resC, 
            LevelData<DATA>& a_phiC, 
            LevelData<DATA>& a_RC,
            LevelData<DATA>& a_phi,
            LevelFluxRegister& a_flux); 
    protected:
        Proto::Stencil<Real> m_flux[DIM];
        Proto::Stencil<Real> m_div[DIM];
        
        Proto::Stencil<Real>        m_avg;
        Proto::InterpStencil<Real>  m_interp;
        Real                        m_dx;
        Real                        m_lambda;
        unsigned int                m_refRatio;
        LevelData<DATA>             m_coarseTemp;
        Proto::InterpStencil<Real>  m_interpBC;
};

template<unsigned int NCOMP, unsigned int NGHOST, class DATA>
void BaseOp<NCOMP, NGHOST, DATA>::define(const DisjointBoxLayout& a_layout,
    Real a_dx, bool a_amr)
{
    m_dx = a_dx;
    if (a_amr)
    {
        m_refRatio = AMR_REFRATIO;
        m_interpBC = Proto::InterpStencil<Real>::Quadratic(m_refRatio);
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
    assert(m_lambda != 0);
    Proto::Point tempGhost;
    if (a_amr)
    {
        tempGhost = this->ghost() + m_interpBC.spanPoint();
    } else {
        tempGhost = this->ghost();
    }
    DisjointBoxLayout coarseLayout;
    coarsen_dbl(coarseLayout, a_layout, m_refRatio);
    m_coarseTemp.define(coarseLayout, NCOMP, tempGhost);
}

template<unsigned int NCOMP, unsigned int NGHOST,class DATA>
void BaseOp<NCOMP, NGHOST, DATA>::apply(patch& a_rhs, patch& a_phi,
        Real a_dx, Real a_scale) const
{
    for (int ii = 0; ii < DIM; ii++)
    {
        patch flx = this->flux(a_phi, ii, a_dx, 1.0);
        if (ii == 0)
        {
            a_rhs |= this->m_div[ii](flx, a_scale/a_dx);
        } else {
            a_rhs += this->m_div[ii](flx, a_scale/a_dx);
        }
    }
}
template<unsigned int NCOMP,  unsigned int NGHOST,class DATA>
void BaseOp<NCOMP, NGHOST, DATA>::apply(LevelData<DATA>& a_rhs,
        LevelData<DATA>& a_phi,
        Real a_dx, Real a_scale) const
{
    a_phi.exchange();
    auto iter = a_phi.dataIterator();
    for (iter.begin(); iter.ok(); ++iter)
    {
        patch phi = a_phi[iter];
        patch rhs = a_rhs[iter];
        apply(rhs, phi, a_dx, a_scale);
    }
}

template<unsigned int NCOMP, unsigned int NGHOST, class DATA>
Real BaseOp<NCOMP, NGHOST, DATA>::residual(LevelData<DATA>& a_res,
        LevelData<DATA>& a_phi,
        const LevelData<DATA>& a_R,
        Real a_dx) const
{
    a_phi.exchange();
    Real resnorm = 0.0;
    auto iter = a_phi.dataIterator();
    for (iter.begin(); iter.ok(); ++iter)
    {
        patch phi = a_phi[iter];
        patch rhs = a_R[iter];
        patch res = a_res[iter];
        apply(res, phi, a_dx, -1.0);
        res += rhs;
        resnorm = std::max(res.absMax(), resnorm);
    }
    return resnorm;
}

template<unsigned int NCOMP, unsigned int NGHOST, class DATA>
Real BaseOp<NCOMP, NGHOST, DATA>::residual(LevelData<DATA>& a_phi,
        const LevelData<DATA>& a_R,
        Real a_dx) const
{
    a_phi.exchange();
    Real resnorm = 0.0;
    auto iter = a_phi.dataIterator();
    for (iter.begin(); iter.ok(); ++iter)
    {
        patch phi = a_phi[iter];
        patch R = a_R[iter];
        patch res(R.box());
        apply(res, phi, a_dx, -1.0);
        res += R;
        resnorm = std::max(res.absMax(), resnorm);
    }
    return resnorm;
}

template<unsigned int NCOMP, unsigned int NGHOST, class DATA>
Real BaseOp<NCOMP, NGHOST, DATA>::residual(LevelData<DATA>& a_res,
        LevelData<DATA>& a_phi,
        const LevelData<DATA>& a_R ) const
{
    return residual(a_res, a_phi, a_R, m_dx);
}

template<unsigned int NCOMP, unsigned int NGHOST, class DATA>
Real BaseOp<NCOMP, NGHOST, DATA>::residual(LevelData<DATA>& a_phi,
        const LevelData<DATA>& a_R ) const
{
    return residual(a_phi, a_R, m_dx);
}

/* UNUSED
template<unsigned int NCOMP, class DATA>
void BaseOp<NCOMP, NGHOST, DATA>::relax(LevelData<DATA>& a_phi,
        LevelData<DATA>& a_res,
        const LevelData<DATA>& a_R,
        int a_n) const
{
    auto iter = a_phi.dataIterator();
    for (int ii = 0; ii < a_n; ii++)
    {
        a_phi.exchange();
        for (iter.begin(); iter.ok(); ++iter)
        {
            patch phi = a_phi[iter];
            patch res = a_res[iter];
            patch R = a_R[iter];
            (*this)(res, phi, m_dx);
            double lambda = m_lambda;
            Proto::forallInPlace(
            [=] PROTO_LAMBDA (var& argPhi, var& argLPhi, const var& argR)
            {
                for (int comp = 0; comp < NCOMP; comp++)
                {
                    argLPhi(comp) = argR(comp) - argLPhi(comp);
                    argPhi(comp) += (argLPhi(comp)*lambda);
                }
            }, phi, res, R);
        } //  end patch iterations
    }// end iterations
}
*/

template<unsigned int NCOMP, unsigned int NGHOST, class DATA>
void BaseOp<NCOMP, NGHOST, DATA>::relax(LevelData<DATA>& a_phi,
        const LevelData<DATA>& a_R,
        int a_n) const
{
    auto layout = a_phi.disjointBoxLayout();
    auto iter = a_phi.dataIterator();
    for (int ii = 0; ii < a_n; ii++)
    {
        a_phi.exchange();
        for (iter.begin(); iter.ok(); ++iter)
        {
            patch phi = a_phi[iter];
            patch R = a_R[iter];
            patch lphi(layout[iter]);
            apply(lphi, phi, m_dx);
            Real lambda = m_lambda;
            Proto::forallInPlace(
            [=] PROTO_LAMBDA (var& argPhi, var& argLPhi, const var& argR)
            {
                for (int comp = 0; comp < NCOMP; comp++)
                {
                    argLPhi(comp) = argR(comp) - argLPhi(comp);
                    argPhi(comp) += (argLPhi(comp)*lambda);
                }
            }, phi, lphi, R);
        } //  end patch iterations
    }// end iterations
}

template<unsigned int NCOMP, unsigned int NGHOST, class DATA>
void BaseOp<NCOMP, NGHOST, DATA>::coarsen(LevelData<DATA>& a_phiC,
        const LevelData<DATA>& a_phi)
{
    auto iter = a_phi.dataIterator();
    for (iter.begin(); iter.ok(); ++iter)
    {
        patch phi = a_phi[iter];
        patch temp = m_coarseTemp[iter];
        temp |= m_avg(phi);
    }
    m_coarseTemp.copyTo(a_phiC);
}

template<unsigned int NCOMP, unsigned int NGHOST, class DATA>
void BaseOp<NCOMP, NGHOST, DATA>::coarseRhs(LevelData<DATA>& a_RC,
        LevelData<DATA>& a_phiC,
        LevelData<DATA>& a_phi,
        const LevelData<DATA>& a_R)
{
    a_phi.exchange();
    a_phiC.exchange();
    // initialize RC with avg of fine residual
    auto fineLayout = a_phi.disjointBoxLayout();
    auto fiter = a_phi.dataIterator();
    for (fiter.begin(); fiter.ok(); ++fiter)
    {
        patch phi = a_phi[fiter];
        patch R = a_R[fiter];
        patch res(fineLayout[fiter]);
        apply(res, phi, m_dx, -1.0);
        res += R;
        patch temp = m_coarseTemp[fiter];
        temp |= m_avg(res);
    }
    m_coarseTemp.copyTo(a_RC);

    // add L(PhiC)
    auto coarseLayout = a_RC.disjointBoxLayout();
    auto citer = a_RC.dataIterator();
    Real cdx = m_dx*m_refRatio;
    for (citer.begin(); citer.ok(); ++citer)
    {
        patch R = a_RC[citer];
        patch phiC = a_phiC[citer];
        patch LPhiC(coarseLayout[citer]);
        apply(LPhiC, phiC, cdx);
        R += LPhiC;
    }
}

template<unsigned int NCOMP, unsigned int NGHOST, class DATA>
void BaseOp<NCOMP, NGHOST, DATA>::fineCorrection(LevelData<DATA>& a_phi,
        LevelData<DATA>& a_phiC,
        const LevelData<DATA>& a_phiC0)
{
    // make phiC0 into a correction
    auto citer = a_phiC.dataIterator();
    for (citer.begin(); citer.ok(); ++citer)
    {
        patch phiC = a_phiC[citer];
        patch phiC0 = a_phiC0[citer];
        phiC0 -= phiC;
        phiC0 *= (-1.0);
    }
    a_phiC0.copyTo(m_coarseTemp);

    auto fineLayout = a_phi.disjointBoxLayout();
    auto fiter = a_phi.dataIterator();
    for (fiter.begin(); fiter.ok(); ++fiter)
    {
        patch phi = a_phi[fiter];
        patch temp = m_coarseTemp[fiter];
        Proto::Box b = Proto::Box(fineLayout[fiter]).coarsen(m_refRatio);
        phi += m_interp(temp, b);
    }
}

template<unsigned int NCOMP, unsigned int NGHOST, class DATA>
void BaseOp<NCOMP, NGHOST, DATA>::interpBoundary(LevelData<DATA>& a_phi,
        const LevelData<DATA>& a_phiC)
{
    a_phiC.copyTo(m_coarseTemp);
    interpBoundary(a_phi, m_coarseTemp, m_interpBC);
}

template<unsigned int NCOMP, unsigned int NGHOST, class DATA>
void BaseOp<NCOMP, NGHOST, DATA>::interpBoundary(LevelData<DATA>& a_phi,
        LevelData<DATA>& a_temp,
        const Proto::InterpStencil<Real>& a_interp)
{
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
        patch coarseData = a_temp[iter];
        patch fineData = a_phi[iter];

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
                fineData |= a_interp(coarseData, intersect);
            }
        }
    }
}

template<unsigned int NCOMP, unsigned int NGHOST, class DATA>
void BaseOp<NCOMP, NGHOST, DATA>::coarseRhs(LevelData<DATA>& a_RC,
        LevelData<DATA>& a_phiC,
        LevelData<DATA>& a_phi,
        const LevelData<DATA>& a_R,
        const LevelData<DATA>& a_rhoC,
        LevelFluxRegister& a_flux)
{
    a_phi.exchange();
    a_phiC.exchange();
    
    auto coarseLayout   = a_phiC.disjointBoxLayout();
    auto fineLayout     = a_phi.disjointBoxLayout();
    Proto::Box coarseDomain = coarseLayout.physDomain().domainBox();
    Proto::Box fineDomain   = fineLayout.physDomain().domainBox();
    auto citer = a_phiC.dataIterator();
    auto fiter = a_phi.dataIterator();
    
    Real cdx = m_dx*m_refRatio;
    Proto::Box invalidRegion = fineDomain.coarsen(m_refRatio);

    a_rhoC.copyTo(a_RC);

    reflux(a_RC, a_phiC, a_phi, a_flux, -1.0);

    for (fiter.begin(); fiter.ok(); ++fiter)
    {
        patch phi = a_phi[fiter];
        patch R = a_R[fiter];
        patch temp = m_coarseTemp[fiter];
        patch res(fineLayout[fiter]);
        apply(res, phi, m_dx, -1.0);
        res += R;
        temp |= m_avg(res);
    }
    m_coarseTemp.copyTo(a_RC);

    for (citer.begin(); citer.ok(); ++citer)
    {
        Proto::Box patchBox = coarseLayout[citer];
        Proto::Box intersect = patchBox & invalidRegion;

        patch phiC = a_phiC[citer];
        patch RC = a_RC[citer];
        patch temp(intersect);
        apply(temp, phiC, cdx);
        RC += temp;
    }
}

template<unsigned int NCOMP, unsigned int NGHOST, class DATA>
void BaseOp<NCOMP, NGHOST, DATA>::reflux(LevelData<DATA>& a_RC,
        LevelData<DATA>& a_phiC,
        LevelData<DATA>& a_phi,
        LevelFluxRegister& a_flux,
        Real a_scale)
{
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
        patch phi = a_phi[fiter];
        for (int dir = 0; dir < DIM; dir++)
        {
            auto fineFlux = FArrayBox(this->flux(phi, dir, m_dx),
                    Proto::Point::Basis(dir));
            a_flux.incrementFine(fineFlux, 1.0, fiter(),
                    comps, comps, dir);
        }
    }
    for (citer.begin(); citer.ok(); ++citer)
    {
        Proto::Box patchBox = coarseLayout[citer];
        patch phiC = a_phiC[citer];
        for (int dir = 0; dir < DIM; dir++)
        {
            auto coarseFlux = FArrayBox(this->flux(phiC, dir, cdx),
                    Proto::Point::Basis(dir));
            a_flux.incrementCoarse(coarseFlux, 1.0, citer(),
                    comps, comps, dir);
        }
    }
    Real scale = a_scale/cdx;
    a_flux.reflux(a_RC, scale);
}

template<unsigned int NCOMP, unsigned int NGHOST, class DATA>
void BaseOp<NCOMP, NGHOST, DATA>::refluxUndivided(LevelData<DATA>& a_RC,
        LevelData<DATA>& a_phiC,
        LevelData<DATA>& a_phi,
        LevelFluxRegister& a_flux,
        Real a_scale)
{
    a_flux.setToZero();
    a_phi.exchange();
    a_phiC.exchange();

    Interval comps(0, NCOMP-1);
   
    auto fineLayout     = a_phi.disjointBoxLayout();
    auto coarseLayout   = a_phiC.disjointBoxLayout();
    auto fiter = a_phi.dataIterator();
    auto citer = a_phiC.dataIterator();
    Proto::Box invalidRegion = fineLayout.physDomain().domainBox();
    invalidRegion = invalidRegion.coarsen(m_refRatio);
    for (fiter.begin(); fiter.ok(); ++fiter)
    {
        Proto::Box patchBox = fineLayout[fiter];
        patch phi = a_phi[fiter];
        for (int dir = 0; dir < DIM; dir++)
        {
            auto fineFlux = FArrayBox(this->flux(phi, dir, 1.0),
                    Proto::Point::Basis(dir));
            a_flux.incrementFine(fineFlux, 1.0, fiter(),
                    comps, comps, dir);
        }
    }
    for (citer.begin(); citer.ok(); ++citer)
    {
        Proto::Box patchBox = coarseLayout[citer];
        patch phiC = a_phiC[citer];
        for (int dir = 0; dir < DIM; dir++)
        {
            auto coarseFlux = FArrayBox(this->flux(phiC, dir, 1.0),
                    Proto::Point::Basis(dir));
            a_flux.incrementCoarse(coarseFlux, 1.0, citer(),
                    comps, comps, dir);
        }
    }
    a_flux.reflux(a_RC, a_scale);

}
/* beginning of a reflux function eliding LevelFluxRegister
template<unsigned int NCOMP, class DATA>
void BaseOp<NCOMP, NGHOST, DATA>::reflux(LevelData<DATA>& a_RC, 
    LevelData<DATA>& a_phi,
    LevelData<DATA>& a_phiC)
{
    auto coarseLayout = a_RC.disjointBoxLayout();
    auto fineLayout = a_phi.disjointBoxLayout();
    auto citer = a_RC.dataIterator();
    auto fiter = a_phi.dataIterator();
    Real cdx = m_dx*m_refRatio;

    Proto::Box invalid = fineLayout.PhysDomain().domainBox();
    invalid = invalid.coarsen(m_refRatio);
    Proto::Box coarseBoundary = invalid.grow(this->ghost());

    LevelData<DATA> FluxDiff[DIM];
    for (int ii = 0; ii < DIM; ii++)
    {
        FluxDiff[ii].define(m_coarseTemp.disjointBoxLayout(), NCOMP, Proto::Point::Ones());
    }

    a_phiC.copyTo(m_coarseTemp);
    
    for (int dir = 0; dir < DIM; dir++)
    {
        for (fiter.begin(); fiter.ok(); ++fiter)
        {
            patch phiC = m_coarseTemp[fiter];
            patch phi = a_phi[fiter];
            patch fluxDiff = FluxDiff[dir][fiter];
            fluxDiff.setVal(0);
            fluxDiff -= this->flux(phiC, dir, cdx);
            fluxDiff += m_avg(this->flux(phi, dir, m_dx));
        }
    }

    for (int ii = 0; ii < DIM; dir++)
    {
        for (fiter.begin(); fiter.ok(); ++fiter)
        {
            patch fluxDiff = FluxDiff[dir][fiter];
            //need to get data back into coarse structure. Not necessarily straightforward
        }
    }
}
*/
template<unsigned int NCOMP, unsigned int NGHOST, class DATA>
void BaseOp<NCOMP, NGHOST, DATA>::coarseResidual(LevelData<DATA>& a_resC, 
        LevelData<DATA>& a_phiC, 
        LevelData<DATA>& a_RC,
        LevelData<DATA>& a_phi,
        LevelFluxRegister& a_flux)
{
    auto citer = a_phiC.dataIterator();
    auto coarseLayout = a_phiC.disjointBoxLayout();
    auto fiter = a_phi.dataIterator();
    auto fineLayout = a_phi.disjointBoxLayout();
    Real cdx = m_dx*m_refRatio;
    Proto::Box invalidRegion = fineLayout.physDomain().domainBox();
    invalidRegion = invalidRegion.coarsen(m_refRatio);

    a_flux.setToZero();
    a_phiC.exchange();
    a_phi.exchange();

    for (citer.begin(); citer.ok(); ++citer)
    {
        patch phiC = a_phiC[citer];
        patch resC = a_resC[citer];
        apply(resC, phiC, cdx, -1.0);
    }
    
    reflux(a_resC, a_phiC, a_phi, a_flux, -1.0);

    for (citer.reset(); citer.ok(); ++citer)
    {
        patch resC = a_resC[citer];
        patch RC = a_RC[citer];
        resC += RC;
        Proto::Box b = Proto::Box(coarseLayout[citer]) & invalidRegion;
        resC.setVal(0,b);
    }
}







