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

        inline BaseOp(){}
        inline void define(const DisjointBoxLayout& a_layout,
                Real a_dx, bool a_amr = true);
        
        inline void buildInterpStencils(int a_maxShift, Proto::Box a_footprint, int a_order);
        virtual inline void init() = 0;

        static constexpr int numcomps() {return NCOMP;} 
        static constexpr Proto::Point ghost() {return Proto::Point::Ones(NGHOST);}
        inline Proto::Point interpGhost() {return this->m_interpBC[0].spanPoint();}
        inline Real dx() const {return m_dx;}
        
        virtual inline patch flux(const patch& a_phi,
                int a_dir, Real a_dx, Real a_scale = 1.0) const = 0;
        virtual inline void apply(patch& a_rhs, patch& a_phi,
                Real a_dx, Real a_scale = 1.0) const;
        virtual inline void apply(LevelData<DATA>& a_rhs,
                LevelData<DATA>& a_phi,
                Real a_dx, Real a_scale = 1.0) const;
        virtual inline void plusApply(LevelData<DATA>& a_rhs,
                LevelData<DATA>& a_phi,
                Real a_dx, Real a_scale = 1.0) const;
        
        virtual inline void residual(patch& a_res,
                patch& a_phi, const patch& a_R, const Real& a_dx) const;
        virtual inline Real residual(LevelData<DATA>& a_res,
                LevelData<DATA>& a_phi,
                const LevelData<DATA>& a_R,
                const Real a_dx) const;
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
                const LevelData<DATA>& a_phiC0,
                int a_n = 0); 

        // AMR functions
        virtual inline void interpBoundary(LevelData<DATA>& a_phi,
                const LevelData<DATA>& a_phiC);
        virtual inline void interpBoundary(LevelData<DATA>& a_phi,
                LevelData<DATA>& a_temp, int a_level);
        virtual inline void reflux(LevelData<DATA>& a_RC,
                LevelData<DATA>& a_phiC,
                LevelData<DATA>& a_phi,
                LevelFluxRegister& a_flux,
                Real a_scale = 1.0);
        virtual inline void coarseResidual(LevelData<DATA>& a_resC,
                const LevelData<DATA>& a_rhoC,
                LevelData<DATA>& a_phiC,
                const LevelData<DATA>& a_R,
                LevelData<DATA>& a_phi,
                LevelFluxRegister& a_flux,
                int a_n = 0);
    protected:
        Proto::Stencil<Real> m_flux[DIM];
        Proto::Stencil<Real> m_div[DIM];
        
        Proto::Stencil<Real>        m_avg;
        Proto::InterpStencil<Real>  m_interp;
        Real                        m_dx;
        Real                        m_lambda;
        unsigned int                m_refRatio;
        LevelData<DATA>             m_coarseTemp;
        std::vector<Proto::InterpStencil<Real>>  m_interpBC;
};

template<unsigned int NCOMP, unsigned int NGHOST, class DATA>
void BaseOp<NCOMP, NGHOST, DATA>::define(const DisjointBoxLayout& a_layout,
    Real a_dx, bool a_amr)
{
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
}

template<unsigned int NCOMP, unsigned int NGHOST,class DATA>
void BaseOp<NCOMP, NGHOST, DATA>::buildInterpStencils(int a_maxShift, Proto::Box a_footprint, int a_order)
{
    int r = MG_REFRATIO;
    while (r <= AMR_REFRATIO)
    {
        m_interpBC.push_back(Proto::InterpStencil<Real>::Build(a_maxShift, a_footprint, a_order, r));
        r *= MG_REFRATIO;
    }
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

template<unsigned int NCOMP,  unsigned int NGHOST,class DATA>
void BaseOp<NCOMP, NGHOST, DATA>::plusApply(LevelData<DATA>& a_rhs,
        LevelData<DATA>& a_phi,
        Real a_dx, Real a_scale) const
{
    a_phi.exchange();
    auto layout = a_phi.disjointBoxLayout();
    auto iter = a_phi.dataIterator();
    for (iter.begin(); iter.ok(); ++iter)
    {
        patch phi = a_phi[iter];
        patch rhs = a_rhs[iter];
        patch temp(layout[iter]);
        apply(temp, phi, a_dx, a_scale);
        rhs += temp;
    }
}

template<unsigned int NCOMP, unsigned int NGHOST, class DATA>
void BaseOp<NCOMP, NGHOST, DATA>::residual(patch& a_res,
        patch& a_phi,
        const patch& a_R,
        const Real& a_dx) const
{
    apply(a_res, a_phi, a_dx, -1.0);
    a_res += a_R;
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
        residual(res, phi, rhs, a_dx);
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
        residual(res, phi, R, a_dx);
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
            patch res(layout[iter]);
           
            residual(res, phi, R, m_dx); 
            res *= m_lambda;
            phi += res;
        }
    }
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
    
    auto fineLayout = a_phi.disjointBoxLayout();
    auto fiter = a_phi.dataIterator();
    for (fiter.begin(); fiter.ok(); ++fiter)
    {
        patch phi = a_phi[fiter];
        patch R = a_R[fiter];
        patch temp = m_coarseTemp[fiter];
        patch res(fineLayout[fiter]);
        
        residual(res, phi, R, m_dx);
        temp |= m_avg(res);
    }
    m_coarseTemp.copyTo(a_RC);
    plusApply(a_RC, a_phiC, m_dx*m_refRatio);
}

template<unsigned int NCOMP, unsigned int NGHOST, class DATA>
void BaseOp<NCOMP, NGHOST, DATA>::fineCorrection(LevelData<DATA>& a_phi,
        LevelData<DATA>& a_phiC,
        const LevelData<DATA>& a_phiC0,
        int a_n)
{
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
    interpBoundary(a_phi, m_coarseTemp, m_interpBC.size() - 1);
}

template<unsigned int NCOMP, unsigned int NGHOST, class DATA>
void BaseOp<NCOMP, NGHOST, DATA>::interpBoundary(LevelData<DATA>& a_phi,
        LevelData<DATA>& a_temp, int a_level)
{
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
                fineData |= interpStencil(coarseData, intersect);
            }
        }
    }
}

template<unsigned int NCOMP, unsigned int NGHOST, class DATA>
void BaseOp<NCOMP, NGHOST, DATA>::reflux(LevelData<DATA>& a_RC,
        LevelData<DATA>& a_phiC,
        LevelData<DATA>& a_phi,
        LevelFluxRegister& a_flux,
        Real a_scale)
{
    //Turn this on to produce debugging outputs from reflux
    bool debug = false;
    
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
    // negative sign -> defaults to refluxing with <F_fine> - F_coarse
    // yes, I am quite sure that this is the way LFR is coded -CLG
    Real scale = -a_scale/cdx;
    a_flux.reflux(a_RC, scale);

    if (debug)
    {
        LevelData<DATA> LPhi, LPhiC, Flux, FluxC_X, FluxC_Y, FluxF_X, FluxF_Y;
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
}

template<unsigned int NCOMP, unsigned int NGHOST, class DATA>
void BaseOp<NCOMP, NGHOST, DATA>::coarseResidual(LevelData<DATA>& a_resC,
        const LevelData<DATA>& a_rhoC,
        LevelData<DATA>& a_phiC,
        const LevelData<DATA>& a_R,
        LevelData<DATA>& a_phi,
        LevelFluxRegister& a_flux,
        int a_n)
{
    Real cdx = m_dx*m_refRatio;
    coarsen(a_phiC, a_phi);
    interpBoundary(a_phi, a_phiC);
    a_phiC.exchange();
    a_phi.exchange();
    
    auto coarseLayout = a_resC.disjointBoxLayout();
    auto citer = a_resC.dataIterator();
    for (citer.begin(); citer.ok(); ++citer)
    {
        patch resC = a_resC[citer];
        patch rhoC = a_rhoC[citer];
        patch phiC = a_phiC[citer];
        apply(resC, phiC, cdx, -1.0);
        resC += rhoC; 
    }
    reflux(a_resC, a_phiC, a_phi, a_flux,1.0);
    
    auto fineLayout = a_phi.disjointBoxLayout();
    auto fiter = a_phi.dataIterator();
    for (fiter.begin(); fiter.ok(); ++fiter)
    {
        patch temp = m_coarseTemp[fiter];
        patch phi = a_phi[fiter];
        patch R = a_R[fiter];
        patch res(fineLayout[fiter]);
        temp.setVal(0);

        residual(res, phi, R, m_dx);
        temp |= m_avg(res);
    }
    m_coarseTemp.copyTo(a_resC);
}
