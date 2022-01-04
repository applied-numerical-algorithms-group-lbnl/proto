#pragma once
#include "Proto.H"

using namespace Proto;

template<typename T, MemType MEM = MEMTYPE_DEFAULT>
class BoxOp_Laplace : public BoxOp<T, 1, 0, MEM>
{
    public:

    // These functions are not optional
    inline static constexpr Point ghost() { return Point::Ones(2);}
    inline static constexpr Point auxGhost() { return Point::Zeros();}
    
    inline BoxOp_Laplace() : BoxOp<T, 1, 0, MEM>() {}
    inline BoxOp_Laplace(double a_dx) : BoxOp<T, 1, 0, MEM>(a_dx) {}

    inline double spectralRadius() const {return (2.0*DIM) / pow(this->dx(), 2); }
    inline void flux(
            BoxData<T>& a_flux,
            const BoxData<T>& a_state,
            int a_dir) const
    {
        auto GRAD = Stencil<T>::DiffCellToFace(a_dir, Side::Lo);
        a_flux |= GRAD(a_state, 1.0 / this->dx()); 
    }
};
