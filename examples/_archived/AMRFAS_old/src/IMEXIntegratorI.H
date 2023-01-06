
template<class IMEXOP, class IMEXDATA, int NSTAGES, int NDENSE>
void IMEXIntegrator<IMEXOP, IMEXDATA, NSTAGES, NDENSE>::init(
    IMEXDATA& a_state,
    double a_dt,
    const std::vector<double>& a_dx,
    bool a_denseOutput)
{
    // define all the coefficients if not done yet. 
    if (!m_defined)
    {
        define();
        m_defined = true;
    }
   
    m_dt = a_dt;
    for (int ii = 0; ii < DIM; ii++)
    {
        m_dx[ii] = a_dx[ii];
    }

    m_rhs.define(a_state);
    m_kE.define(a_state);
    m_kI.define(a_state);
    m_phi.resize(NSTAGES);
    for (int istage = 0; istage < NSTAGES; istage++)
    {
        m_phi[istage] = std::make_shared<IMEXDATA>(a_state);
    }
    if (a_denseOutput)
    {
        m_denseCoefs.resize(NDENSE);
        for (int icoef = 0; icoef < NDENSE; icoef++)
        {
            m_denseCoefs[icoef] = std::make_shared<IMEXDATA>(a_state);
        }
    }
    m_denseOutput = a_denseOutput;
    m_hasDenseOutput = false;
    m_init = true;
}
template<class IMEXOP, class IMEXDATA, int NSTAGES, int NDENSE>
void IMEXIntegrator<IMEXOP, IMEXDATA, NSTAGES, NDENSE>::advance(double a_time, IMEXDATA& a_state)
{
    PROTO_ASSERT(m_init, "Tried to call IMEXIntegrator::advance before it was initialized.");
    if (m_denseOutput)
    {
        m_hasDenseOutput = false;
        for (int icoef = 0; icoef < NDENSE; icoef++)
        {
            m_denseCoefs[icoef]->clear();
        }
    }
    m_op.preProcess(a_state);
    for (int istage = 0; istage < NSTAGES; istage++)
    {
        a_state.copyTo(*(m_phi[istage]));
    }
    
    for (int istage = 0; istage < NSTAGES; istage++)
    {
        double t = a_time + m_c[istage]*m_dt;
        if (istage > 0)
        {
            m_phi[istage]->copyTo(m_rhs);
            m_op.solve(*m_phi[istage], m_rhs, t, istage);
        }
        
        m_op.explicitOp(m_kE, *m_phi[istage], t, istage);
        a_state.increment(m_kE, m_dt*m_b[istage]);
        for (int k = istage+1; k < NSTAGES; k++)
        {
            m_phi[k]->increment(m_kE, m_dt*m_aE[k][istage]);
        }
        
        m_op.implicitOp(m_kI, *m_phi[istage], t, istage);
        a_state.increment(m_kI, m_dt*m_b[istage]);
        for (int k = istage+1; k < NSTAGES; k++)
        {
            m_phi[k]->increment(m_kI, m_dt*m_aI[k][istage]);
        }
        
        if (m_denseOutput)
        {
            for (int icoef = 0; icoef < NDENSE; icoef++)
            {
                m_denseCoefs[icoef]->increment(m_kE, m_dt*m_bstar[icoef][istage]);
                m_denseCoefs[icoef]->increment(m_kI, m_dt*m_bstar[icoef][istage]);
            }
        }
        m_hasDenseOutput = m_denseOutput;
    }
    m_op.postProcess(a_state);
}





