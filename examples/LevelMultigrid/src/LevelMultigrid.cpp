#include "LevelMultigrid.H"
using namespace std;

LevelMultigrid::LevelMultigrid(
        const DisjointBoxLayout& a_layout,
        double a_dx,
        int a_level
        )
{
    this->define(a_layout,a_dx,a_level);
}

void LevelMultigrid::define(
        const DisjointBoxLayout& a_layout,
        double a_dx,
        int a_level)
{
    PR_TIMERS("LevelMultigrid::define");
    m_layout = a_layout;

    m_level = a_level;
    m_dx = a_dx;
    m_lambda = m_dx*m_dx/(4*DIM);
    Point boxsize = m_layout.boxSize();
    ProblemDomain pd = m_layout.domain();
    if (m_level > 0)
    {
        Point refRatio = Point::Ones(2);
        PROTO_ASSERT(m_layout.coarsenable(refRatio),
                "LevelMultigrid::define | Error: Layout is not coarsenable by refinment ratio. \
                You may be using too many levels for your problem size.");
        DisjointBoxLayout dblCoarseLocal = m_layout.coarsen(refRatio);
        DisjointBoxLayout dblCoarse;
        ProblemDomain pdCoarse = pd.coarsen(refRatio);
        if (pdCoarse.sizes() % boxsize == Point::Zeros())
        {
            dblCoarse = DisjointBoxLayout(pdCoarse,boxsize);
        }
        else
        {
            dblCoarse =  DisjointBoxLayout(pdCoarse,pdCoarse.sizes());
        }
        m_resc.define(dblCoarse,Point::Zeros());
        m_delta.define(dblCoarse,Point::Ones());
        m_localCoarse.define(dblCoarseLocal, Point::Zeros());
        m_coarsePtr =
            shared_ptr<LevelMultigrid>(new LevelMultigrid(dblCoarse,2*m_dx,m_level-1));
    }
}

void LevelMultigrid::coarseResidual(
        LevelBoxData<double >& a_resc,
        LevelBoxData<double >& a_phi,
        LevelBoxData<double >& a_rhs)
{
    PR_TIMERS("LevelMultigrid::coarseResidual");
    a_phi.exchange();
    double hsqinv = 1./(m_dx*m_dx);
    auto L = Stencil<double>::Laplacian();
    auto A = Stencil<double>::AvgDown(2);
    for (auto dit : a_phi)
    {
        BoxData<double>& phi = a_phi[dit];
        BoxData<double>& rhs = a_rhs[dit];
        BoxData<double>& rescLocal = m_localCoarse[dit];
        BoxData<double> res(a_phi.layout()[dit]);
        res.setVal(0.);
        res += rhs;
        res += L(phi,-hsqinv);
        rescLocal |= A(res);
    }
    m_localCoarse.copyTo(a_resc);
}

double LevelMultigrid::resnorm(
        LevelBoxData<double >& a_phi,
        LevelBoxData<double >& a_rhs)
{
    PR_TIMERS("LevelMultigrid::resnorm");
    a_phi.exchange();
    double hsqinv = 1./(m_dx*m_dx);

    m_rxn.reset();
    auto L = Stencil<double>::Laplacian();
    for (auto dit : a_phi)
    {
        BoxData<double>& phi = a_phi[dit];
        BoxData<double>& rhs = a_rhs[dit];
        BoxData<double> res(a_phi.layout()[dit]);

        res.setVal(0.);
        res -= rhs;
        res += L(phi,hsqinv);
        res.reduce(m_rxn);
    }
    return m_rxn.fetch();
}

void LevelMultigrid::pointRelax(
        LevelBoxData<double >& a_phi,
        LevelBoxData<double >& a_rhs,
        int a_numIter)
{
    PR_TIMERS("LevelMultigrid::pointRelax");
    double weight = 1.0/(4*DIM);
    Stencil<double> D = (-m_lambda)*Shift(Point::Zeros());
    auto L = Stencil<double>::Laplacian();
    m_rxn.reset();
    for (int iter = 0; iter < a_numIter; iter++)
    {
        a_phi.exchange();
        for (auto dit : a_phi)
        {
            BoxData<double>& phi = a_phi[dit];
            BoxData<double>& rhs = a_rhs[dit];
            BoxData<double> temp = L(phi, weight);
            temp += D(rhs); 
            if (iter == a_numIter - 1) { temp.reduce(m_rxn);}
            phi += temp;
        }
    }
}

void LevelMultigrid::fineInterp(
        LevelBoxData<double >& a_phi,
        LevelBoxData<double >& a_delta)
{
    PR_TIMERS("LevelMultigrid::fineInterp");
    a_delta.copyTo(m_localCoarse);

    for (auto dit : a_phi)
    {
        BoxData<double>& phi = a_phi[dit];
        BoxData<double>& delta = m_localCoarse[dit];

        Box K(Point::Zeros(),Point::Ones());
        phi += m_fineInterp(delta);
    }
}

void LevelMultigrid::vCycle(
        LevelBoxData<double >& a_phi,
        LevelBoxData<double >& a_rhs)
{
    PR_TIMERS("LevelMultigrid::vCycle");
    if (m_level > 0)
    {
        pointRelax(a_phi,a_rhs,m_preRelax);
        coarseResidual(m_resc,a_phi,a_rhs);
        m_delta.setToZero();
        m_coarsePtr->vCycle(m_delta,m_resc);
        fineInterp(a_phi,m_delta);
        pointRelax(a_phi,a_rhs,m_postRelax);
    }
    else
    {
        pointRelax(a_phi,a_rhs,m_bottomRelax);
    }
}
