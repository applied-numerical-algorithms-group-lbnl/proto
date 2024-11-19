#pragma once
#include "Proto.H"

using namespace Proto;

template<typename T,
    template<MemType> typename MAP,
    MemType MEM = MEMTYPE_DEFAULT>
class BoxOp_MBLaplace : public MBBoxOp<T, MAP, 1, 0, MEM>
{
    public:

    // inherit constructors
    using MBBoxOp<T, MAP, 1, 0, MEM>::MBBoxOp;
    
    // How many ghost cells does the operator need from the state variables
    inline static Point ghost() { return Point::Ones(3);}
    
    // How many ghost cells does the operator need from the auxiliary variables
    inline static Point auxGhost() { return Point::Zeros();}
    
    // What is the intended order of accuracy of the operator
    inline static constexpr int order() { return 4; }

    // The approximate spectral radius. Often needed if the operator is to be inverted using 
    // an iterative solver
    inline double spectralRadius() const {return (2.0*DIM) / pow(this->dxMin(), 2); }
#if 0
    inline void
    operator()(
            BoxData<T,1>& a_output,
            const BoxData<T,1>& a_state,
            T                a_scale = 1.0) const
    {
        Array<BoxData<T,1>, DIM> fluxes;
        for (int dir = 0; dir < DIM; dir++)
        {
            fluxes[dir].define(a_output.box().grow(dir, Side::Hi, 1));
        }
        (*this)(a_output, fluxes, a_state, a_scale);
    }
#endif
    inline void
    apply(
            BoxData<T,1>& a_output,
            Array<BoxData<T,1>, DIM>& a_fluxes,
            const BoxData<T,1>& a_state) const override
    {
        for (int dd = 0; dd < DIM; dd++)
        {
            flux(a_fluxes[dd], a_state, dd);
        }
        //TODO: BCs are not passing tests yet. 
        //this->boundaryConditions(a_fluxes, a_state);
        this->divergence(a_output, a_fluxes);
    }

    inline void applyBC (
            Array<BoxData<T,1>, DIM>& a_fluxes,
            const BoxData<T,1>& a_state,
            Face a_face) const override
    {
        if (!m_useBC) { return; }
        auto _dx = this->dx();
        Box domainBox = this->layout().getBlock(this->block()).domain().box();
        BoundaryCondition::Dirichlet(a_fluxes, a_state, 0.0, a_face, domainBox, _dx);
    }

    inline void
    flux(
            BoxData<T, 1>& a_flux,
            const BoxData<T, 1>& a_state,
            int a_dir) const
    {
        auto faceGrad =  Operator::_faceGradxPhi(
                a_state, a_state, m_NTMatrix[a_dir], m_NTMatrix[a_dir],
                m_Jf[a_dir], m_Jf[a_dir], a_dir);
        auto tmp = Operator::_faceMatrixProductAB(
                faceGrad,m_NT[a_dir],faceGrad,m_NT[a_dir],a_dir);
        tmp.copyTo(a_flux);
        a_flux /= this->dx()[a_dir];
    }

    // Cache metric terms and Jacobian
    inline void init()
    {
#if DIM == 2
        m_NT.define(this->x().box().extrude(Point::Ones(),-1));
#elif DIM == 3
        m_NT.define(this->x().box().grow(-2).extrude(Point::Ones(),-1));
#else
        MayDay<void>::Error("BoxOp_MBLaplace not implemented for DIM != 2,3");
#endif
        m_NTMatrix.clear();
        m_Jf.clear();
        for (int dir = 0; dir < DIM; dir++)
        {
            auto C2F = Stencil<double>::CellToFace(dir);
            // In 2D, NT does not require any ghost cells from X
            // In 3D, NT requires 2 layers of ghost cells from X
            auto NT_dir = Operator::cofactor(this->x(), dir, this->dx()); 
            NT_dir.copyTo(m_NT[dir]);
            m_Jf.push_back(C2F(this->jacobian()));
        }
        for (int dir = 0; dir < DIM; dir++)
        {
            //NT should have (conservative) 2 layers of ghost cells
            m_NTMatrix.push_back(Operator::cofactorMatrix(m_NT, dir));
        }
        m_useBC = true;
    }

    inline void toggleBC(bool flag) { m_useBC = flag; }

    private:
    bool m_useBC;
    FluxBoxData<T,DIM> m_NT;
    std::vector<BoxData<T,DIM,MEM,DIM>> m_NTMatrix;
    std::vector<BoxData<T,1,MEM>> m_Jf;
};
