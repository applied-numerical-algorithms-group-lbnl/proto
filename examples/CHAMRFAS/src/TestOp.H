#pragma once

#include "Proto.H"
#include "LevelData.H"
#include "AMRUtils.H"
#include "DisjointBoxLayout.H"
#include "LevelFluxRegister.H"
#include "Proto_WriteBoxData.H"
#include <cmath>
#define NUMCOMPS 1

//using namespace Proto;

typedef Proto::BoxData<Real, NUMCOMPS> BD;
typedef Proto::Var<Real, NUMCOMPS> V;
typedef BaseFab<Real> FAB;

template<class DATA>
class TestOp
{
    public:
        
        typedef Proto::BoxData<Real, NUMCOMPS> patchData;
        
        inline TestOp(){}
        inline TestOp(const DisjointBoxLayout& a_layout, Real a_dx){define(a_layout, a_dx);};
        inline void define(const DisjointBoxLayout& a_layout, Real a_dx, bool a_amr = true);

        static constexpr int numcomps(){return NUMCOMPS;}

        inline Proto::Point ghost() const {return Proto::Point::Ones();}  
        inline Real dx() const {return m_dx;};
        
        inline Proto::BoxData<Real> flux(const BD& a_phi,int dir, Real a_dx = 0, Real a_scale = 1.0) const;
        inline Proto::BoxData<Real> div(const BD& a_flux, int dir, Real a_dx = 0, Real a_scale = 1.0) const;
        inline void operator()(BD& a_rhs, const BD& a_phi, Real a_dx = 0, Real a_scale = 1.0) const;
        inline void operator()(LevelData<DATA>& a_rhs, const LevelData<DATA>& a_phi, Real a_dx = 0, Real a_scale = 1.0) const;
       
        Real residual(LevelData<DATA>& a_res,
                LevelData<DATA>& a_phi,
                const LevelData<DATA>& a_R,
                Real a_dx = 0) const;
        Real residual(LevelData<DATA>& a_phi,
                const LevelData<DATA>& a_R,
                Real a_dx = 0) const;

        // Multigrid stuff
        void relax(LevelData<DATA>& a_phi,
                LevelData<DATA>& a_res,
                const LevelData<DATA>& a_R,
                int a_n) const;
        void relax(LevelData<DATA>& a_phi,
                const LevelData<DATA>& a_R,
                int a_n) const;
        void coarsen(LevelData<DATA>& a_phiC,
                LevelData<DATA>& a_phi);
        void coarseRhs(LevelData<DATA>& a_RC,
                LevelData<DATA>& a_phiC,
                LevelData<DATA>& a_phi,
                const LevelData<DATA>& a_R);
        void fineCorrection(LevelData<DATA>& a_phi,
                LevelData<DATA>& a_phiC,
                const LevelData<DATA>& a_phiC0);

        // AMR stuff
        void updateInterp(LevelData<DATA>& a_phi, const LevelData<DATA>& a_phiC);
        void interpBoundary(LevelData<DATA>& a_phi, 
                const LevelData<DATA>& a_phiC);
        void interpBoundary(LevelData<DATA>& a_phi, 
                LevelData<DATA>& a_temp,
                const Proto::InterpStencil<Real>& a_interp);
        void computeRhs(LevelData<DATA>& a_RC, 
                LevelData<DATA>& a_phiC,
                LevelData<DATA>& a_phi,
                const LevelData<DATA>& a_R,
                const LevelData<DATA>& a_rhoC,
                LevelFluxRegister& a_flux);
        void reflux(LevelData<DATA>& a_RC, 
                LevelData<DATA>& a_phiC,
                LevelData<DATA>& a_phi,
                LevelFluxRegister& a_flux);

    private:
        Proto::Stencil<Real> m_avg;            //averaging stencil
        Proto::Stencil<Real> m_flux[DIM];      
        Proto::Stencil<Real> m_div[DIM];
        Proto::InterpStencil<Real> m_interp;   //interpolation stencil
        Real m_dx;                      
        Real m_lambda;
        unsigned int m_refRatio;

        LevelData<DATA> m_coarseTemp;   //temp for coarsening / refining in multigrid

        Proto::InterpStencil<Real> m_interpBC;   //AMR BC interpolation stencil

};

template<class DATA>
void TestOp<DATA>::define(const DisjointBoxLayout& a_layout, Real a_dx, bool a_amr)
{
    CH_TIME("TestOp::define");
    Proto::Point tempGhost;
    if (a_amr) 
    {
        tempGhost = Proto::Point::Ones(2);
        m_refRatio = AMR_REFRATIO;
    }
    else 
    {
        tempGhost = Proto::Point::Ones();
        m_refRatio = MG_REFRATIO;
    }
#if DIM == 2
    //TODO: Are these still right?
    m_lambda = -0.5*(a_dx*a_dx)*(6.0/20.0);
#elif DIM == 3
    //TODO: Are these still right?
    m_lambda = -0.5*(a_dx*a_dx)*(15.0/64.0);
#endif          //undefined otherwise for the time being
    m_dx = a_dx;
    m_avg = Proto::Stencil<Real>::AvgDown(m_refRatio);
    // operator stencil is a 2*DIM + 1 Point Laplacian
    for (int ii = 0; ii < DIM; ii++)
    {
        m_flux[ii] = 1.0*Proto::Shift::Zeros() - 1.0*Proto::Shift::Basis(ii,-1);
        m_div[ii] = 1.0*Proto::Shift::Basis(ii,1) - 1.0*Proto::Shift::Zeros();
    }
    m_interp = Proto::InterpStencil<Real>::PiecewiseConstant(Proto::Point::Ones(m_refRatio));
    m_interpBC = Proto::InterpStencil<Real>::Quadratic(m_refRatio);

    DisjointBoxLayout coarseLayout;
    coarsen_dbl(coarseLayout, a_layout, m_refRatio);
    m_coarseTemp.define(coarseLayout, NUMCOMPS, tempGhost);
}

template<class DATA>
Proto::BoxData<Real> TestOp<DATA>::flux(const BD& a_phi,int dir, Real a_dx, Real a_scale) const
{
    Real dx = m_dx;
    if (a_dx > 0) {dx = a_dx;}
    Proto::BoxData<Real> ret = m_flux[dir](a_phi, a_scale/dx);
    return ret;
} 

// Compute divergence (of a flux)
template<class DATA>
Proto::BoxData<Real> TestOp<DATA>::div(const BD& a_flux, int dir, Real a_dx, Real a_scale) const
{
    Real dx = m_dx;
    if (a_dx > 0) {dx = a_dx;}
    
    Proto::BoxData<Real> ret;
    ret = m_div[dir](a_flux,a_scale/dx);
    return ret;
} 

template<class DATA>
void TestOp<DATA>::operator()(BD& a_rhs, const BD& a_phi, Real a_dx, Real a_scale) const
{
    Real dx = m_dx;
    if (a_dx > 0) {dx = a_dx;}
    
    for (int ii = 0; ii < DIM; ii++)
    {
        BD flux = m_flux[ii](a_phi, 1.0/dx);
        if (ii == 0)
        {
            a_rhs |= m_div[ii](flux, a_scale/dx);
        } else {
            a_rhs += m_div[ii](flux, a_scale/dx);
        }
    }
}

template<class DATA>
void TestOp<DATA>::operator()(LevelData<DATA>& a_rhs, const LevelData<DATA>& a_phi, Real a_dx, Real a_scale) const
{
    Real dx = m_dx;
    if (a_dx > 0) {dx = a_dx;}
 
    auto iter = a_phi.dataIterator();
    for (iter.begin(); iter.ok(); ++iter)
    {
        BD phi = a_phi[iter];
        BD rhs = a_rhs[iter];
        (*this)(rhs, phi, dx, a_scale);
    }
}

template<class DATA>
Real TestOp<DATA>::residual(LevelData<DATA>& a_res,
        LevelData<DATA>& a_phi,
        const LevelData<DATA>& a_rhs,
        Real a_dx) const
{
    CH_TIME("TestOp::residual");
    a_phi.exchange();

    Real dx = m_dx;
    if (a_dx > 0) {dx = a_dx;}

    Real resnorm = 0.0;    
    auto iter = a_phi.dataIterator();
    for (iter.begin(); iter.ok(); ++iter)
    {
        BD phi = a_phi[iter];
        BD rhs = a_rhs[iter];
        BD res = a_res[iter];
        (*this)(res, phi, dx, -1.0);
        res += rhs;
        resnorm = std::max(res.absMax(), resnorm);
    }
    return resnorm;
}

template<class DATA>
Real TestOp<DATA>::residual(LevelData<DATA>& a_phi,
        const LevelData<DATA>& a_R,
        Real a_dx) const
{
    CH_TIME("TestOp::residual");
    a_phi.exchange();
    
    Real dx = m_dx;
    if (a_dx > 0) {dx = a_dx;}

    Real resnorm = 0.0;
    auto iter = a_phi.getBoxes().dataIterator();
    for (iter.begin(); iter.ok(); ++iter)
    {
        BD phi = a_phi[iter()];
        BD R = a_R[iter()];
        BD res(R.box());
        (*this)(res, phi, dx, -1.0);
        res += R;
        resnorm = std::max(res.absMax(), resnorm);
    }
    return resnorm;
}

template<class DATA>
void TestOp<DATA>::relax(LevelData<DATA>& a_phi,
        LevelData<DATA>& a_res,
        const LevelData<DATA>& a_R, int a_n) const
{
    CH_TIMERS("TestOp::relax");
    
    auto iter = a_phi.dataIterator();
    for (int ii = 0; ii < a_n; ii++)
    {
        a_phi.exchange();
        for (iter.begin(); iter.ok(); ++iter)
        {
            BD phi = a_phi[iter()];
            BD res = a_res[iter()];
            BD R = a_R[iter()];
            (*this)(res, phi, m_dx);
            double lambda = m_lambda;
            Proto::forallInPlace(
                [lambda] PROTO_LAMBDA (V& a_phi, V& a_lphi, const V& a_R)
                { 
                    for (int comp = 0; comp < NUMCOMPS; comp++)
                    {
                        a_lphi(comp) = a_R(comp) - a_lphi(comp);
                        a_phi(comp) += (a_lphi(comp)*lambda);
                    }
                }, phi, res, R);
        }
    }
}

template<class DATA>
void TestOp<DATA>::relax(LevelData<DATA>& a_phi,
        const LevelData<DATA>& a_R, int a_n) const
{
    CH_TIMERS("TestOp::relax");

    auto layout = a_phi.disjointBoxLayout();
    auto iter = a_phi.dataIterator();
    for (int ii = 0; ii < a_n; ii++)
    {
        a_phi.exchange();
        for (iter.begin(); iter.ok(); ++iter)
        {
            BD phi = a_phi[iter()];
            BD R = a_R[iter()];
            BD lphi(layout[iter]);
            (*this)(lphi, phi, m_dx);
            double lambda = m_lambda;
            Proto::forallInPlace(
                [lambda] PROTO_LAMBDA (V& a_phi, V& a_lphi, const V& a_R)
                { 
                    for (int comp = 0; comp < NUMCOMPS; comp++)
                    {
                        a_lphi(comp) = a_R(comp) - a_lphi(comp);
                        a_phi(comp) += (a_lphi(comp)*lambda);
                    }
                }, phi, lphi, R);
        }
    }
}

template<class DATA>
void TestOp<DATA>::coarsen(LevelData<DATA>& a_phiC,
        LevelData<DATA>& a_phi)
{
    CH_TIME("TestOp::coarsen");
    //a_phi.exchange();
    auto iter = a_phi.getBoxes().dataIterator();
    for (iter.begin(); iter.ok(); ++iter)
    {
        BD phi = a_phi[iter()];
        BD temp = m_coarseTemp[iter()];
        temp |= m_avg(phi);
    }
    m_coarseTemp.copyTo(a_phiC);
}

// Multigrid Version
template<class DATA>
void TestOp<DATA>::coarseRhs(LevelData<DATA>& a_RC,
        LevelData<DATA>& a_phiC,
        LevelData<DATA>& a_phi,
        const LevelData<DATA>& a_R)
{
    CH_TIME("TestOp::coarseRhs");
    a_phi.exchange(); 
    auto fiter = a_phi.getBoxes().dataIterator();
    for (fiter.begin(); fiter.ok(); ++fiter)
    {
        BD phi = a_phi[fiter()];
        BD R = a_R[fiter()];
        BD res(R.box());
        (*this)(res, phi, m_dx, -1.0);
        res += R;
        BD temp = m_coarseTemp[fiter()];
        temp |= m_avg(res);
    }
    m_coarseTemp.copyTo(a_RC);  
    a_phiC.exchange(); 
    auto citer = a_RC.getBoxes().dataIterator();

    for (citer.begin(); citer.ok(); ++citer)
    {
        BD R = a_RC[citer()];
        BD phiC = a_phiC[citer()];
        double cdx = m_dx*m_refRatio;
        BD LPhiC(R.box());
        (*this)(LPhiC, phiC, cdx);
        R += LPhiC;
    }
}

template<class DATA>
void TestOp<DATA>::fineCorrection(LevelData<DATA>& a_phi,
        LevelData<DATA>& a_phiC,
        const LevelData<DATA>& a_phiC0)
{
    CH_TIME("TestOp::fineCorrection");
    
    // Compute coarse correction
    auto citer = a_phiC.dataIterator();
    for (citer.begin(); citer.ok(); ++citer)
    {
        BD phiC = a_phiC[citer()];
        BD phiC0 = a_phiC0[citer()];
        
        phiC0 -= phiC; 
        phiC0 *= (-1); //PhiC0 is now a correction
    }
    a_phiC0.copyTo(m_coarseTemp);
    
    // Interpolate correction
    auto fiter = a_phi.dataIterator();
    for (fiter.begin(); fiter.ok(); ++fiter)
    {
        BD phi = a_phi[fiter()];
        BD temp = m_coarseTemp[fiter()];
        phi += m_interp(temp);
    }
}

//========================================================================
// AMR Stuff
//========================================================================
template<typename DATA>
void TestOp<DATA>::updateInterp(LevelData<DATA>& a_phi, const LevelData<DATA>& a_phiC)
{
    // get references to layouts and iterators
    auto iter = a_phi.dataIterator(); //for a_phi and m_coarseTemp
    auto fineLayout = a_phi.disjointBoxLayout();
    auto tempLayout = m_coarseTemp.disjointBoxLayout();
    auto coarseLayout = a_phiC.disjointBoxLayout();
   
    // move coarse data over to temp array
    a_phiC.copyTo(m_coarseTemp);
    
    // update interior data with averaged fine values
    for (iter.begin(); iter.ok(); ++iter)
    {
        BD coarse = m_coarseTemp[iter()];
        BD fine = a_phi[iter()];
        coarse |= m_avg(fine,tempLayout[iter()]);
    }
    
    // fix up the internal boundaries
    m_coarseTemp.exchange();
    interpBoundary(a_phi, m_coarseTemp, m_interpBC);
}

template<typename DATA>
void TestOp<DATA>::interpBoundary(LevelData<DATA>& a_phi,
                                  const LevelData<DATA>& a_phiC)
{
    a_phiC.copyTo(m_coarseTemp);
    interpBoundary(a_phi, m_coarseTemp, m_interpBC);
}


template<typename DATA>
void TestOp<DATA>::interpBoundary(LevelData<DATA>& a_phi,
                                  LevelData<DATA>& a_temp,
                                  const Proto::InterpStencil<Real>& a_interp)
{
    auto iter = a_phi.dataIterator();
    auto fineLayout = a_phi.disjointBoxLayout();
    Proto::Box fineDomain = fineLayout.physDomain().domainBox();
    
    // Compute bitRatio and bitBox
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
        BD coarseData = a_temp[iter];
        BD fineData = a_phi[iter];

        Proto::Box coarsePatch = coarseData.box();
        Proto::Point bitPoint = Proto::Box(fineLayout[iter]).low()/bitRatio;
        for (auto niter = U.begin(); niter != U.end(); ++niter)
        {
            if (*niter == Proto::Point::Zeros()) {continue;}
            Proto::Point bitNeighbor = bitPoint + *niter;
            if (!bitBox.contains(bitNeighbor))
            {
                Proto::Box coarseNeighborBox(bitNeighbor, bitNeighbor);
                coarseNeighborBox = coarseNeighborBox.refine(bitRatio/m_refRatio);
                Proto::Box intersect = coarseNeighborBox & coarsePatch;
                fineData |= a_interp(coarseData, intersect);
            }
        }
    }
}

template<typename DATA>
void TestOp<DATA>::computeRhs(LevelData<DATA>& a_RC, 
        LevelData<DATA>& a_phiC,
        LevelData<DATA>& a_phi,
        const LevelData<DATA>& a_R,
        const LevelData<DATA>& a_rhoC,
        LevelFluxRegister& a_flux)
{
    Real cdx = m_dx*m_refRatio;
    
    // get handles for legibility
    a_phi.exchange();
    a_phiC.exchange();
    const DisjointBoxLayout& coarseLayout = a_phiC.disjointBoxLayout();
    const DisjointBoxLayout& fineLayout   = a_phi.disjointBoxLayout();
    Proto::Box coarseDomain = coarseLayout.physDomain().domainBox();
    Proto::Box fineDomain   = fineLayout.physDomain().domainBox();
    auto citer = a_phiC.dataIterator();
    auto fiter = a_phi.dataIterator();
    
    // determine the invalid region of the coarse level
    Proto::Box invalidRegion = fineDomain.coarsen(m_refRatio);

    // initialize with rho
    a_rhoC.copyTo(a_RC);
    
    // compute reflux divergence
    reflux(a_RC, a_phiC, a_phi, a_flux); //reflux the data
    
    // replace invalid region with averaged down residual
    
    for (fiter.begin(); fiter.ok(); ++fiter)
    {
        BD phi = a_phi[fiter];
        BD R = a_R[fiter];
        BD temp = m_coarseTemp[fiter];
        BD res(R.box());
        (*this)(res, phi, m_dx, -1.0);
        res += R;
        temp |= m_avg(res);
    }
    m_coarseTemp.copyTo(a_RC); // *should* overwrite garbage from refluxing
    
    // Add L(phi_coarse)
    for (citer.begin(); citer.ok(); ++citer)
    {
        Proto::Box patchBox = coarseLayout[citer];
        Proto::Box intersect = patchBox & invalidRegion;
        BD phiC = a_phiC[citer];
        BD RC = a_RC[citer];
        BD temp(intersect);
        (*this)(temp, phiC, cdx);
        RC += temp;
    }
}

// Just do refluxing
template<typename DATA>
void TestOp<DATA>::reflux(LevelData<DATA>& a_RC, 
        LevelData<DATA>& a_phiC,
        LevelData<DATA>& a_phi,
        LevelFluxRegister& a_flux)
{
    a_flux.setToZero();
    a_phi.exchange();
    a_phiC.exchange();
    const DisjointBoxLayout& coarseLayout = a_phiC.disjointBoxLayout();
    const DisjointBoxLayout& fineLayout = a_phi.disjointBoxLayout();

    Interval comps(0,0);

    auto citer = a_phiC.dataIterator();
    auto fiter = a_phi.dataIterator();
    
    for (fiter.begin(); fiter.ok(); ++fiter)
    {
        Proto::Box patchBox = fineLayout[fiter];
        const BD phi = a_phi[fiter];
        for (int dir = 0; dir < DIM; dir++)
        {
            auto fluxF = FArrayBox(this->flux(phi, dir, m_dx), Proto::Point::Basis(dir)); //compute flux
            a_flux.incrementFine(fluxF, 1.0, fiter(), comps, comps, dir); //increment flux register
        }
    }
    
    for (citer.begin(); citer.ok(); ++citer)
    {
        Proto::Box patchBox = coarseLayout[citer];
        const BD phiC = a_phiC[citer];
        for (int dir = 0; dir < DIM; dir++)
        {
            auto fluxC = FArrayBox(this->flux(phiC, dir, m_dx*m_refRatio), Proto::Point::Basis(dir)); //compute flux
          
            a_flux.incrementCoarse(fluxC, 1.0, citer(), comps, comps, dir); //increment flux register
        }
    }
   
    Real scale = -1.0/(m_refRatio*m_dx); //coarse dx
    a_flux.reflux(a_RC, scale);
}
