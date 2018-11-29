#pragma once

#include "Proto.H"
#include "LevelData.H"
#include "AMRUtils.H"
#include "DisjointBoxLayout.H"
#include <cmath>
#define NUMCOMPS 1

//using namespace Proto;

int fileNum = 0;
typedef Proto::BoxData<Real, NUMCOMPS> BD;
typedef Proto::Var<Real, NUMCOMPS> V;
typedef BaseFab<Real> FAB;

template<class DATA>
class TestOp
{
    public:


        inline TestOp(){}
        inline TestOp(const DisjointBoxLayout& a_layout, Real a_dx){define(a_layout, a_dx);};
        inline void define(const DisjointBoxLayout& a_layout, Real a_dx);

        static constexpr int numcomps(){return NUMCOMPS;}

        inline Proto::Point ghost() const {return m_L.ghost();}  
        inline Real dx() const {return m_dx;};

        // Multigrid stuff
        Real residual(LevelData<DATA>& a_res,
                LevelData<DATA>& a_phi,
                const LevelData<DATA>& a_R) const;
        Real residual(LevelData<DATA>& a_phi,
                const LevelData<DATA>& a_R) const;
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
                const LevelData<DATA>& a_res);
        void coarseRhs(LevelData<DATA>& a_RC,
                LevelData<DATA>& a_phiC,
                const LevelData<DATA>& a_phi,
                const LevelData<DATA>& a_R);
        void fineCorrection(LevelData<DATA>& a_phi,
                LevelData<DATA>& a_phiC,
                const LevelData<DATA>& a_phiC0);

        // AMR stuff
        void interpBoundary(LevelData<DATA>& a_phi, const LevelData<DATA>& a_phiC);

        // Initialization
        inline static void initialCondition(LevelData<DATA>& a_phi, Real dx);
        inline static void forcing(LevelData<DATA>& a_R, Real dx);
        inline static void solution(LevelData<DATA>& a_soln, Real dx);

    private:
        Proto::Stencil<Real> m_L;              //actual operator
        Proto::Stencil<Real> m_avg;            //averaging stencil
        Proto::InterpStencil<Real> m_interp;   //interpolation stencil
        Real m_dx;                      
        Real m_lambda;

        LevelData<DATA> m_coarseTemp;   //temp for coarsening / refining in multigrid

        Proto::InterpStencil<Real> m_interpBC;   //AMR BC interpolation stencil

};

    template<class DATA>
void TestOp<DATA>::define(const DisjointBoxLayout& a_layout, Real a_dx)
{
    CH_TIME("TestOp::define");
#if DIM == 2
    m_L = Proto::Stencil<Real>::Laplacian_9();
    m_lambda = -0.5*(a_dx*a_dx)*(6.0/20.0);
#elif DIM == 3
    m_L = Proto::Stencil<Real>::Laplacian_27();
    m_lambda = -0.5*(a_dx*a_dx)*(15.0/64.0);
#endif          //undefined otherwise for the time being
    m_dx = a_dx;
    m_avg = Proto::Stencil<Real>::AvgDown(MG_REFRATIO);
    m_interp = Proto::InterpStencil<Real>::PiecewiseConstant(Proto::Point::Ones(MG_REFRATIO));
    m_interpBC = Proto::InterpStencil<Real>::Quadratic(AMR_REFRATIO);
    DisjointBoxLayout coarseLayout;
    coarsen_dbl(coarseLayout, a_layout, MG_REFRATIO);
    // extra ghost data for interpolating boundary conditions;
    m_coarseTemp.define(coarseLayout, NUMCOMPS, Proto::Point::Ones(2));
}

template<class DATA>
Real TestOp<DATA>::residual(LevelData<DATA>& a_res,
        LevelData<DATA>& a_phi,
        const LevelData<DATA>& a_R) const
{
    CH_TIME("TestOp::residual");
    a_phi.exchange();
    auto iter = a_phi.getBoxes().dataIterator();
    Real resnorm = 0;
    for (iter.begin(); iter.ok(); ++iter)
    {
        BD phi = a_phi[iter()];
        BD res = a_res[iter()];
        BD R = a_R[iter()];
        res |= m_L(phi,1.0/(m_dx*m_dx));
        Proto::forallInPlace(
                [](V& a_lphi, const V& a_R)
                {
                for (int comp = 0; comp < NUMCOMPS; comp++)
                {
                a_lphi(comp) = a_R(comp) - a_lphi(comp);
                }
                }, res, R);
        resnorm = std::max(res.absMax(), resnorm);
    }
    return resnorm;
}

template<class DATA>
Real TestOp<DATA>::residual(LevelData<DATA>& a_phi,
        const LevelData<DATA>& a_R) const
{
    CH_TIME("TestOp::residual");
    a_phi.exchange();
    auto iter = a_phi.getBoxes().dataIterator();
    Real resnorm = 0;
    for (iter.begin(); iter.ok(); ++iter)
    {
        BD phi = a_phi[iter()];
        BD R = a_R[iter()];
        BD res = m_L(phi,1.0/(m_dx*m_dx));
        Proto::forallInPlace(
                [](V& a_lphi, const V& a_R)
                {
                for (int comp = 0; comp < NUMCOMPS; comp++)
                {
                a_lphi(comp) = a_R(comp) - a_lphi(comp);
                }
                }, res, R);
        resnorm = std::max(res.absMax(), resnorm);
    }
    return resnorm;
}

template<class DATA>
void TestOp<DATA>::relax(LevelData<DATA>& a_phi,
        LevelData<DATA>& a_res,
        const LevelData<DATA>& a_R, int a_n) const
{
    //TODO: make sure LD inputs all have same DBL
    CH_TIMERS("TestOp::relax");
    CH_TIMER("TestOp::relax::operator",t1);
    CH_TIMER("TestOp::relax::forall",t2);

    auto iter = a_phi.getBoxes().dataIterator();
    for (int ii = 0; ii < a_n; ii++)
    {
        a_phi.exchange();
        for (iter.begin(); iter.ok(); ++iter)
        {
            BD phi = a_phi[iter()];
            BD res = a_res[iter()];
            BD R = a_R[iter()];
            CH_START(t1);
            res |= m_L(phi, 1.0/(m_dx*m_dx));
            CH_STOP(t1);
            double lambda = m_lambda;
            CH_START(t2);
            Proto::forallInPlace(
                [lambda](V& a_phi, V& a_lphi, const V& a_R)
                { 
                    for (int comp = 0; comp < NUMCOMPS; comp++)
                    {
                        a_lphi(comp) = a_R(comp) - a_lphi(comp);
                        a_phi(comp) += (a_lphi(comp)*lambda);
                    }
                }, phi, res, R);
            CH_STOP(t2);
        }
    }
}

template<class DATA>
void TestOp<DATA>::relax(LevelData<DATA>& a_phi,
        const LevelData<DATA>& a_R, int a_n) const
{
    //TODO: make sure LD inputs all have same DBL
    CH_TIMERS("TestOp::relax");
    CH_TIMER("TestOp::relax::operator",t1);
    CH_TIMER("TestOp::relax::forall",t2);

    auto iter = a_phi.getBoxes().dataIterator();
    for (int ii = 0; ii < a_n; ii++)
    {
        a_phi.exchange();
        for (iter.begin(); iter.ok(); ++iter)
        {
            BD phi = a_phi[iter()];
            BD R = a_R[iter()];
            CH_START(t1);
            BD res = m_L(phi, 1.0/(m_dx*m_dx));
            //std::cout << res.absMax() << std::endl;
            CH_STOP(t1);
            double lambda = m_lambda;
            CH_START(t2);
            Proto::forallInPlace(
                [lambda](V& a_phi, V& a_lphi, const V& a_R)
                { 
                    for (int comp = 0; comp < NUMCOMPS; comp++)
                    {
                        a_lphi(comp) = a_R(comp) - a_lphi(comp);
                        a_phi(comp) += (a_lphi(comp)*lambda);
                    }
                }, phi, res, R);
            CH_STOP(t2);
        }
    }
}

    template<class DATA>
void TestOp<DATA>::coarsen(LevelData<DATA>& a_phiC,
        LevelData<DATA>& a_phi)
{
    //TODO: make sure LD inputs all have same DBL
    CH_TIME("TestOp::coarsen");
    a_phi.exchange();
    auto iter = a_phi.getBoxes().dataIterator();
    for (iter.begin(); iter.ok(); ++iter)
    {
        BD phi = a_phi[iter()];
        BD temp = m_coarseTemp[iter()];
        //Box b = m_coarseTemp.getBoxes()[iter()];
        temp |= m_avg(phi);
    }
    m_coarseTemp.copyTo(a_phiC);
}


template<class DATA>
void TestOp<DATA>::coarseRhs(LevelData<DATA>& a_RC,
        LevelData<DATA>& a_phiC,
        const LevelData<DATA>& a_res)
{
    CH_TIME("TestOp::coarseRhs");
    auto fiter = a_res.getBoxes().dataIterator();
    for (fiter.begin(); fiter.ok(); ++fiter)
    {
        BD res = a_res[fiter()];
        BD temp = m_coarseTemp[fiter()];
        //Box b = m_coarseTemp.getBoxes()[fiter()];
        temp |= m_avg(res);
    }
    m_coarseTemp.copyTo(a_RC);  
    a_phiC.exchange(); 
    auto citer = a_RC.getBoxes().dataIterator();

    for (citer.begin(); citer.ok(); ++citer)
    {
        BD R = a_RC[citer()];
        BD phiC = a_phiC[citer()];
        //Box b = a_RC.getBoxes()[citer()];
        double dx = m_dx*MG_REFRATIO;
        R += m_L(phiC,1.0/(dx*dx));
    }
}

template<class DATA>
void TestOp<DATA>::coarseRhs(LevelData<DATA>& a_RC,
        LevelData<DATA>& a_phiC,
        const LevelData<DATA>& a_phi,
        const LevelData<DATA>& a_R)
{
    CH_TIME("TestOp::coarseRhs");
    auto fiter = a_phi.getBoxes().dataIterator();
    for (fiter.begin(); fiter.ok(); ++fiter)
    {
        BD phi = a_phi[fiter()];
        BD R = a_R[fiter()];
        BD res = m_L(phi,1.0/(m_dx*m_dx));
        Proto::forallInPlace(
            [](V& a_lphi, const V& a_R)
            {
                for (int comp = 0; comp < NUMCOMPS; comp++)
                {
                    a_lphi(comp) = a_R(comp) - a_lphi(comp);
                }
            }, res, R);
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
        //Box b = a_RC.getBoxes()[citer()];
        double dx = m_dx*MG_REFRATIO;
        R += m_L(phiC,1.0/(dx*dx));
    }
}

    template<class DATA>
void TestOp<DATA>::fineCorrection(LevelData<DATA>& a_phi,
        LevelData<DATA>& a_phiC,
        const LevelData<DATA>& a_phiC0)
{
    CH_TIME("TestOp::fineCorrection");
    auto citer = a_phiC.dataIterator();
    for (citer.begin(); citer.ok(); ++citer)
    {
        BD phiC = a_phiC[citer()];
        BD phiC0 = a_phiC0[citer()];
        phiC0 -= phiC; //phiC0 is now a correction
        phiC0 *= (-1);
    }
    a_phiC0.copyTo(m_coarseTemp);

    auto fiter = a_phi.getBoxes().dataIterator();
    for (fiter.begin(); fiter.ok(); ++fiter)
    {
        BD phi = a_phi[fiter()];
        BD temp = m_coarseTemp[fiter()];
        //Box b = m_coarseTemp.getBoxes()[fiter()];
        phi += m_interp(temp);
    }
}

    template<class DATA>
void TestOp<DATA>::initialCondition(LevelData<DATA>& a_phi, Real dx)
{
    CH_TIME("TestOp::initialCondition");
    auto iter = a_phi.dataIterator();
    for (iter.begin(); iter.ok(); ++iter)
    {
        BD phi = a_phi[iter()];
        Proto::forallInPlace_p([dx](Proto::Point p, V& phi)
                {
                double x = p[0]*dx + 0.5*dx + .25*M_PI;
                phi(0) = -cos(x) + 0.1*cos(x);
                phi(0) = 0.;
                },phi);
    }
}

    template<class DATA>
void TestOp<DATA>::forcing(LevelData<DATA>& a_R, Real dx)
{
    CH_TIME("TestOp::forcing");
    auto iter = a_R.dataIterator();
    for (iter.begin(); iter.ok(); ++iter)
    {
        BD R = a_R[iter()];
        Proto::forallInPlace_p([dx](Proto::Point p, V& R)
                {
                double x = p[0]*dx + 0.5*dx;
                R(0) = cos(x);
                },R);
    }
}

    template<class DATA>
void TestOp<DATA>::solution(LevelData<DATA>& a_soln, Real dx)
{
    CH_TIME("TestOp::solution");
    auto iter = a_soln.dataIterator();
    for (iter.begin(); iter.ok(); ++iter)
    {
        BD soln = a_soln[iter()];
        Proto::forallInPlace_p([dx](Proto::Point p, V& soln)
                {
                double x = p[0]*dx + 0.5*dx;
                soln(0) = -cos(x);
                },soln);
    }
}

//========================================================================
// AMR Stuff
//========================================================================
    template<typename DATA>
void TestOp<DATA>::interpBoundary(LevelData<DATA>& a_phi, const LevelData<DATA>& a_phiC)
{
    // move coarse data over to temp array
    a_phiC.copyTo(m_coarseTemp);
    // average down the fine data where it exists
    auto iter = a_phi.dataIterator();
    auto fineLayout = a_phi.getBoxes();
    auto coarseLayout = m_coarseTemp.getBoxes();
    for (iter.begin(); iter.ok(); ++iter)
    {
        BD coarse = m_coarseTemp[iter()];
        BD fine = a_phi[iter()];
        coarse |= m_avg(fine,coarseLayout[iter()]);
    }
    // fix up the internal boundaries
    m_coarseTemp.exchange();
    // compute a bit box
    Box domainBox = coarseLayout.physDomain().domainBox();
    int r[DIM];
    for (int ii = 0; ii < DIM; ii++)
    {
        //pretty sure this logic only works with power of 2 refinement -clg
        r[ii] = std::min((int)domainBox.size(ii), MAXBOXSIZE);
    }
    Proto::Point bitRatio(r);
    Proto::Box bitbox = domainBox.coarsen(bitRatio);

    // do interpolation
    Proto::Box U = Proto::Box::Cube(3).shift(-Proto::Point::Ones());
    for (iter.begin(); iter.ok(); ++iter)
    {
        Proto::Point bitPoint = Proto::Box(coarseLayout[iter]).low()/bitRatio;
        for (auto niter = U.begin(); niter != U.end(); ++niter)
        {
            if (*niter == Proto::Point::Zeros()){continue;}
            Proto::Point neighbor = bitPoint + *niter;
            if (!bitbox.contains(neighbor)) //this is an external neighbor point
            {
                BD coarseData = m_coarseTemp[iter];
                BD fineData = a_phi[iter];
                Proto::Box neighborBox(neighbor, neighbor);
                neighborBox = neighborBox.refine(bitRatio);
                Proto::Box intersect = neighborBox & coarseData.box();
                fineData |= m_interpBC(coarseData, intersect);
            }
        }
    }
}
