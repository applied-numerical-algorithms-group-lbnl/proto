
template<class OP>
void Multigrid<OP>::define(const DisjointBoxLayout& a_layout,
        double     a_dx,
        int        a_level,
        bool       a_amr,
        int        a_amrLevel)  
{
    PR_TIME("Multigrid::define");
    bool debug = false;
    m_level = a_level;
    
    m_op.define(a_layout, a_dx, false);
    if (a_amr && (a_amrLevel > 0))
    { 
        unsigned int num_mg = log(AMR_REFRATIO)/log(MG_REFRATIO);
        //std::cout << "number of MG levels: " << num_mg << std::endl;
        unsigned int interpRatio = AMR_REFRATIO/ipow(MG_REFRATIO, num_mg - m_level - 1); 
        if (debug) {
            std::cout << "Building Multigrid Level: " << m_level << std::endl;
            std::cout << "\tAMR = true. AMR Level = " << a_amrLevel << std::endl;
            std::cout << "\tNumber of Multigrid sublevels = " << num_mg << std::endl;
            std::cout << "\tRatio between this level and next coarsest AMR level: " << interpRatio << std::endl;
        }
        
        DisjointBoxLayout coarseAMRLayout;
        coarsen_dbl(coarseAMRLayout, a_layout, interpRatio);
        m_phiCAMR.define(coarseAMRLayout, OP::numcomps(), m_op.ghost() + m_op.interpGhost());
    }
    if (a_level > 0)
    {
        DisjointBoxLayout coarseLayout;
        coarsenGrid(coarseLayout, a_layout, MG_REFRATIO);

        m_phiC0.define(coarseLayout, OP::numcomps(), OP::ghost());
        m_phiC.define(coarseLayout, OP::numcomps(), OP::ghost());
        m_RC.define(coarseLayout, OP::numcomps(), Proto::Point::Zeros());

        m_coarser = make_shared<Multigrid<OP>>(coarseLayout, a_dx*MG_REFRATIO, a_level, a_amr, a_amrLevel);
    }
    m_defined = true;    
}

template<class OP>
void Multigrid<OP>::vcycle(LD& a_phi, const LD& a_R)
{
    PR_TIME("Multigrid::vcycle (vanilla Multigrid)");
    bool debug = false;
    if (debug) {
        writeLevel(a_phi, "MG_BOTTOM_phi_L%i_0.hdf5", m_level);
        std::cout << "\tIntegral of Phi: MG_BOTTOM Level " << m_level << ", MGVcycle start: " << integrate(a_phi, m_op.dx()) << std::endl;
    }
    if (m_level == 0)
    {
        m_op.relax(a_phi, a_R, BOTTOM_RELAX);
    } else {
        m_op.relax(a_phi, a_R, PRE_RELAX);
        if (debug) {
            writeLevel(a_phi, "MG_BOTTOM_phi_L%i_1.hdf5", m_level);
            std::cout << "\tIntegral of Phi: MG_BOTTOM Level " << m_level << ", Down relax: " << integrate(a_phi, m_op.dx()) << std::endl;
        }
        m_op.coarsen(m_phiC, a_phi);
        m_phiC.copyTo(m_phiC0);
        m_op.coarseRhs(m_RC, m_phiC, a_phi, a_R);
        m_coarser->vcycle(m_phiC, m_RC);
        if (debug) {
            writeLevel(a_phi, "MG_BOTTOM_phi_L%i_2.hdf5", m_level);
            std::cout << "\tIntegral of Phi: MG_BOTTOM Level " << m_level << ", Recursive call: " << integrate(a_phi, m_op.dx()) << std::endl;
        }
        m_op.fineCorrection(a_phi, m_phiC, m_phiC0);
        m_op.relax(a_phi, a_R, POST_RELAX);
        if (debug) {
            writeLevel(a_phi, "MG_BOTTOM_phi_L%i_3.hdf5", m_level);
            std::cout << "\tIntegral of Phi: MG_BOTTOM Level " << m_level << ", Up relax: " << integrate(a_phi, m_op.dx()) << std::endl;
        }
    }
}

template<class OP>
void Multigrid<OP>::vcycle(LD& a_phi, LD& a_phiC, const LD& a_R)
{
    PR_TIME("Multigrid::vcycle (AMR Multigrid)");
    bool debug = false;
    a_phiC.copyTo(m_phiCAMR);
    if (debug) {
        writeLevel(m_phiCAMR, "MG_phiCAMR.hdf5");
        writeLevel(a_phi, "MG_phi_L%i_0.hdf5", m_level);
        std::cout << "\tIntegral of Phi: MGLevel " << m_level << ", MGVcycle start: " << integrate(a_phi, m_op.dx()) << std::endl;
    }
    m_op.interpBoundary(a_phi, m_phiCAMR, m_level);
    if (debug) {
        writeLevel(a_phi, "MG_phi_L%i_1.hdf5", m_level);
        std::cout << "\tIntegral of Phi: MGLevel " << m_level << ", After Boundary Interp: " << integrate(a_phi, m_op.dx()) << std::endl;
    }
    if (m_level == 0)
    {
        m_op.relax(a_phi, a_R, BOTTOM_RELAX);
        if (debug) {
            writeLevel(a_phi, "MG_phi_L%i_2.hdf5", m_level);
            std::cout << "\tIntegral of Phi: MGLevel " << m_level << ", After bottom relax: " << integrate(a_phi, m_op.dx()) << std::endl;
        }
    } else {
        m_op.relax(a_phi, a_R, PRE_RELAX);
        if (debug) {
            writeLevel(a_phi, "MG_phi_L%i_2.hdf5", m_level);
            std::cout << "\tIntegral of Phi: MGLevel " << m_level << ", After pre relax: " << integrate(a_phi, m_op.dx()) << std::endl;
        }
        m_op.coarsen(m_phiC, a_phi);
        a_phiC.copyTo(m_coarser->m_phiCAMR);
        m_coarser->m_op.interpBoundary(m_phiC, m_coarser->m_phiCAMR, m_level-1);
        m_phiC.copyTo(m_phiC0);
        m_op.coarseRhs(m_RC, m_phiC, a_phi, a_R);
        
        m_coarser->vcycle(m_phiC, a_phiC, m_RC);
        
        m_op.fineCorrection(a_phi, m_phiC, m_phiC0);
        if (debug) {
            writeLevel(a_phi, "MG_phi_L%i_3.hdf5", m_level);
            std::cout << "\tIntegral of Phi: MGLevel " << m_level << ", After fine correction: " << integrate(a_phi, m_op.dx()) << std::endl;
        }
        m_op.relax(a_phi, a_R, POST_RELAX);
        if (debug) {
            writeLevel(a_phi, "MG_phi_L%i_4.hdf5", m_level);
            std::cout << "\tIntegral of Phi: MGLevel " << m_level << ", After post relax: " << integrate(a_phi, m_op.dx()) << std::endl;
        }
    }
}
