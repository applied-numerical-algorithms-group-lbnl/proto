template<class AMR_OP, class DATA>
inline void 
AMRFAS<AMR_OP, DATA>::define(const std::vector<DisjointBoxLayout>& a_grids,
        double a_dx,
        int    a_level,
        int    a_bottomDepth)
{
    m_level = a_level;
    m_op.define(a_grids[a_level], a_dx);
    if (a_level > 0)
    {
        m_mg.define(a_grids[a_level], a_dx, AMR_REFRATIO/MG_REFRATIO - 1, true);
        m_phiC0.define(a_grids[a_level-1], AMR_OP::numcomps(), m_mg.op().ghost());
        m_RC.define(a_grids[a_level-1], AMR_OP::numcomps(), Proto::Point::Zeros());
        m_coarser = make_shared<AMRFAS<AMR_OP, DATA>>(a_grids, a_dx*AMR_REFRATIO, a_level - 1, a_bottomDepth);
        m_reflux.define(a_grids[a_level], a_grids[a_level-1],
                a_grids[a_level].physDomain(), AMR_REFRATIO, 
                AMR_OP::numcomps());
    }
    else
    {
        m_mg.define(a_grids[a_level], a_dx, a_bottomDepth, false);
    }
    m_defined = true;
}

template<class AMR_OP, class DATA>
inline void
AMRFAS<AMR_OP, DATA>::vcycle(std::vector<std::shared_ptr<LevelData<DATA>>> a_phi,
        const std::vector<std::shared_ptr<LevelData<DATA>>> a_rho,
        const std::vector<std::shared_ptr<LevelData<DATA>>> a_res)
{
    this->vcycle(a_phi, a_rho, a_res, *(a_rho.back()));
}

template<class AMR_OP, class DATA>
inline void
AMRFAS<AMR_OP, DATA>::vcycle(std::vector<std::shared_ptr<LevelData<DATA>>> a_phi,
        const std::vector<std::shared_ptr<LevelData<DATA>>> a_rho,
        const std::vector<std::shared_ptr<LevelData<DATA>>> a_res,
        const LevelData<DATA>& a_R)
{
    AMR_LEVEL = m_level;
    // alias values
    LevelData<DATA>& phi = *(a_phi[m_level]);    
    LevelData<DATA>& res = *(a_res[m_level]);
    m_op.residual(res, phi, a_R);
    write(a_phi, "Phi.L%i.0.hdf5", m_level);
    write(a_res, "Res.L%i.0.hdf5", m_level);
    writeLevel(a_R, "R.L%i.0.hdf5", m_level);
    if (m_level == 0)
    {
        auto iter = res.dataIterator();
        m_op.residual(res, phi, a_R);
        
        // no coarser data exists; no need to interp BCs
        m_mg.vcycle(phi,a_R); //bottom solve
        m_op.residual(res, phi, a_R);

        write(a_phi, "Phi.L%i.1.hdf5", m_level);
        write(a_res, "Res.L%i.1.hdf5", m_level);
    }
    else
    {
        LevelData<DATA>& phiC = *(a_phi[m_level-1]);
        LevelData<DATA>& rhoC = *(a_rho[m_level-1]);
        
        m_op.interpBoundary(phi, phiC); 
        m_mg.vcycle(phi, phiC, a_R);
        m_op.residual(res, phi, a_R);
        write(a_res, "Res.L%i.1.hdf5", m_level);
        write(a_phi, "Phi.L%i.1.hdf5", m_level);
        
        m_op.coarsen(phiC, phi);
        phiC.copyTo(m_phiC0);
        write(a_phi, "Phi.L%i.2.hdf5", m_level);
       
        if (m_level > 1) //at least 2 coarser levels exist
        {
            LevelData<DATA>& phiCC = *(a_phi[m_level-2]);
            m_coarser->m_op.updateInterp(phiC, phiCC);
        } 
        
        m_op.computeRhs(m_RC, phiC, phi, a_R, rhoC, m_reflux);
        writeLevel(m_RC, "RC.L%i.0.hdf5", m_level);

        m_coarser->vcycle(a_phi, a_rho, a_res, m_RC);
        m_op.residual(res, phi, a_R);
        write(a_phi, "Phi.L%i.3.hdf5", m_level);
        write(a_res, "Res.L%i.2.hdf5", m_level);
       
        m_op.fineCorrection(phi, phiC, m_phiC0);
        m_op.residual(res, phi, a_R);
        write(a_phi, "Phi.L%i.4.hdf5", m_level);
        write(a_res, "Res.L%i.3.hdf5", m_level);
       
        m_op.interpBoundary(phi, phiC); 
        m_op.residual(res, phi, a_R);
        write(a_phi, "Phi.L%i.5.hdf5", m_level);
        write(a_res, "Res.L%i.4.hdf5", m_level);
        
        m_mg.vcycle(phi, phiC, a_R);
        m_op.residual(res, phi, a_R);
        write(a_phi, "Phi.L%i.6.hdf5", m_level);
        write(a_res, "Res.L%i.5.hdf5", m_level);
    }
}

template<class AMR_OP, class DATA>
inline void
AMRFAS<AMR_OP, DATA>::operator()(std::vector<std::shared_ptr<LevelData<DATA>>> a_rhs,
        std::vector<std::shared_ptr<LevelData<DATA>>> a_phi)
{
    LevelData<DATA>& phi = *(a_phi[m_level]);    
    LevelData<DATA>& rhs = *(a_rhs[m_level]);
    if (m_level > 0)
    {   
        LevelData<DATA>& phiC = *(a_phi[m_level-1]);    
        LevelData<DATA>& rhsC = *(a_rhs[m_level-1]);    
        m_op.coarsen(phiC, phi);
        (*m_coarser)(a_rhs, a_phi);
        m_op.reflux(rhsC, phiC, phi, m_reflux);
        writeLevel(phi, "Phi.op.0.hdf5");
        writeLevel(phiC, "PhiC.op.0.hdf5");
        m_op.interpBoundary(phi, phiC);
        writeLevel(phi, "Phi.op.1.hdf5");
        m_op(rhs, phi);
        m_op.coarsen(rhsC, rhs); //to fix garbage from reflux in invalid regions
    } else {
        m_op(rhs, phi);
    }
}

template<class AMR_OP, class DATA>
template<typename... Args>
inline void
AMRFAS<AMR_OP, DATA>::write(std::vector<std::shared_ptr<LevelData<FArrayBox>>> a_data,
    const char* a_fname, Args... a_params)
{
    char fname[100];
    sprintf(fname, a_fname, a_params...);
    Vector<DisjointBoxLayout> grids(m_level+1);
    Vector<LevelData<FArrayBox>*> data(m_level+1);
    for (int ii = 0; ii < grids.size(); ii++)
    {
        grids[ii] = (*(a_data[ii])).disjointBoxLayout();
        auto& temp = *(a_data[ii]);
        data[ii] = &temp;
    }

    Vector<string> varNames(NUMCOMPS);
    for (int ii = 0; ii < varNames.size(); ii++)
    {
        char name[100];
        sprintf(name, "var_%i", ii);
        std::string s(name);
        varNames[ii] = s;
    }

    Box domain = (*(a_data[0])).disjointBoxLayout().physDomain().domainBox();
    Real dx = m_mg.op().dx()*ipow(AMR_REFRATIO, m_level);
    Real dt = 1.0;
    Real t = 0.0;
    Vector<int> ratios(m_level, AMR_REFRATIO);
    int numLevels = m_level + 1;
    WriteAMRHierarchyHDF5(fname, grids, data, varNames, domain, dx, dt, t, ratios, numLevels);
}
