#pragma once
#include "Proto.H"

using namespace Proto;

template<typename T, MemType MEM = MEMTYPE_DEFAULT>
class BoxOp_TestLaplace : public BoxOp<T, 1, 0, MEM>
{
    public:

    // inherit constructors
    using BoxOp<T, 1, 0, MEM>::BoxOp;
    
    // How many ghost cells does the operator need from the state variables
    inline static Point ghost() { return Point::Ones(2);}
    
    // What is the intended order of accuracy of the operator
    inline static constexpr int order() { return 4; }

    inline void flux(
            BoxData<T, 1>& a_flux,
            const BoxData<T, 1>& a_state,
            int a_dir) const
    {
        auto GRAD = Stencil<T>::DiffCellToFace(a_dir, Side::Lo);
        a_flux |= GRAD(a_state, 1.0 / this->dx()[a_dir]); 
    }

    private:
};
