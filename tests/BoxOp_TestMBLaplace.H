#pragma once
#include "Proto.H"

using namespace Proto;

template<typename T, MemType MEM = MEMTYPE_DEFAULT>
class BoxOp_TestMBLaplace : public BoxOp<T, 1, 0, MEM>
{
    public:

    // inherit constructors
    using BoxOp<T, 1, 0, MEM>::BoxOp;
    
    // How many ghost cells does the operator need from the state variables
    inline static Point ghost() { return Point::Ones(1);}
    
    // How many ghost cells does the operator need from the auxiliary variables
    inline static Point auxGhost() { return Point::Ones(1);}
    
    // What is the intended order of accuracy of the operator
    inline static constexpr int order() { return 2; }

    // The approximate spectral radius. Often needed if the operator is to be inverted using 
    // an iterative solver
    inline double spectralRadius() const {return 1;}
   
    inline void flux(
            BoxData<T, 1>& a_flux,
            const BoxData<T, 1>& a_state,
            int a_dir) const
    {
        auto flux =  Operator::_faceGradxPhi(
                a_state, a_state, m_NTMatrix[a_dir], m_NTMatrix[a_dir],
                m_Jf[a_dir], m_Jf[a_dir], a_dir);
        a_flux = Operator::_faceMatrixProductAB(
                flux,m_NT[a_dir],flux,m_NT[a_dir],a_dir);
    }

    inline void init()
    {
        m_NT.clear();
        m_NTMatrix.clear();
        m_Jf.clear();
        Array<BoxData<T,DIM>, DIM> NT;
        for (int dir = 0; dir < DIM; dir++)
        {
            m_NT.push_back(Operator::cofactor(this->x(),dir));
            m_Jf.push_back(Stencil<double>::CellToFace(dir)(this->jacobian()));
            NT[dir].define(m_NT[dir].box());
            m_NT[dir].copyTo(NT[dir]);
        }
        for (int dir = 0; dir < DIM; dir++)
        {
            m_NTMatrix.push_back(Operator::cofactorMatrix(NT, dir));
        }
    }

    private:

    std::vector<BoxData<T, DIM>> m_NT;
    std::vector<BoxData<T,DIM,MEM,DIM>> m_NTMatrix;
    std::vector<BoxData<T,1,MEM>> m_Jf;
};
