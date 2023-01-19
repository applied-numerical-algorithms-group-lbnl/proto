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
        return Operator::_faceGradxPhi(
                a_state, a_state, m_NTMatrix, m_NTMatrix, m_Jf, m_Jf, a_dir);
    }

    inline void init()
    {
        for (int dir = 0; dir < DIM; dir++)
        {
            m_NT[dir] = Operator::cofactor(this->x(),dir);
            m_Jf[dir] = Stencil<double>::CellToFace(dir)(this->jacobian());
        }
        for (int dir = 0; dir < DIM; dir++)
        {
            m_NTMatrix = Operator::_cofactorMatrix(m_NT, dir); 
        }
    }

    private:

    Array<BoxData<T, DIM>, DIM> m_NT;
    Array<BoxData<T,DIM,MEM,DIM>, DIM> m_NTMatrix;
    Array<BoxData<T,1,MEM>, DIM> m_Jf;
};
