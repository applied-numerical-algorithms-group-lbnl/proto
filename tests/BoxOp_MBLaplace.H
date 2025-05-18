#pragma once
#include "Proto.H"

using namespace Proto;

template<typename T,
    typename MAP,
    MemType MEM = MEMTYPE_DEFAULT>
class BoxOp_MBLaplace : public MBBoxOp<T, MAP, 1, 0, MEM>
{
    public:

    // inherit constructors
    using MBBoxOp<T, MAP, 1, 0, MEM>::MBBoxOp;
    
    // How many ghost cells does the operator need from the state variables
    inline static Point ghost()
    {
        Point diffGhost = Point::Zeros();
        Point matProdGhost = Point::Zeros();
        for (int ii = 0; ii < DIM; ii++)
        {
            diffGhost = absMaxPoint(diffGhost, MappedFaceAverageDiffOp<T>(ii).ghost());
            matProdGhost = absMaxPoint(matProdGhost, FaceAverageMatrixProductOp<T>(ii).ghost());
        }
        return Point::Ones() + diffGhost + matProdGhost;
    }
    
    // How many ghost cells does the operator need from the auxiliary variables
    inline static Point auxGhost() { return Point::Zeros();}
    
    // What is the intended order of accuracy of the operator
    inline static constexpr int order() { return 4; }

    // The approximate spectral radius. Often needed if the operator is to be inverted using 
    // an iterative solver
    inline double spectralRadius() const {return (2.5*DIM) / pow(this->dxMin(), 2); }

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
            BoxData<T, 1>& flux,
            const BoxData<T, 1>& state,
            int dir) const
    {
        auto faceGrad = (*gradOps[dir])(state, state, adjJ_f[dir], adjJ_f[dir], J_f[dir], J_f[dir]);
        (*matProdOps[dir])(flux, faceGrad, N_f[dir], faceGrad, N_f[dir]);
        flux /= this->dx()[dir];
    }

    // Cache metric terms and Jacobian
    inline void init()
    {
        auto areas = this->da();
    
        for (int dir = 0; dir < DIM; dir++)
        {
            gradOps[dir] = std::make_unique<MappedFaceAverageDiffOp<T>>(dir);
            matProdOps[dir] = std::make_unique<FaceAverageMatrixProductOp<T>>(dir);

            auto cellToFaceOp = Stencil<T>::CellToFace(dir);
            FaceAverageAdjugateMatrixOp<T> adjOp(dir);
            FaceAverageSurfaceTransformOp<T,DIM> surfOp(dir);
            J_f[dir] = cellToFaceOp(this->jacobian());
            N_f[dir] = surfOp(this->x());
            N_f[dir] /= areas[dir];
            adjJ_f[dir] = adjOp(this->x());
            adjJ_f[dir] /= areas[dir];
        }
        m_useBC = true;
    }

    inline void toggleBC(bool flag) { m_useBC = flag; }

    private:
    bool m_useBC;
    Array<std::unique_ptr<MappedFaceAverageDiffOp<T>>,DIM> gradOps;
    Array<std::unique_ptr<FaceAverageMatrixProductOp<T>>,DIM> matProdOps;
    Array<BoxData<T,DIM>,DIM> N_f;
    Array<BoxData<T,DIM,MEM,DIM>,DIM> adjJ_f;
    Array<BoxData<T,1,MEM>,DIM> J_f;
};
