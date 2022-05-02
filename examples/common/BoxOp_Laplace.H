#pragma once
#include "Proto.H"

using namespace Proto;

template<typename T, MemType MEM = MEMTYPE_DEFAULT>
class BoxOp_Laplace : public BoxOp<T, 1, 0, MEM>
{
    public:

    // These functions are not optional
    inline static Point ghost() { return Point::Ones(2);}
    inline static Point auxGhost() { return Point::Zeros();}
    inline static constexpr int order() { return 4; }

    inline BoxOp_Laplace() : BoxOp<T, 1, 0, MEM>() {}
    inline BoxOp_Laplace(double a_dx) : BoxOp<T, 1, 0, MEM>(a_dx) {}
    inline BoxOp_Laplace(std::array<double, DIM> a_dx) : BoxOp<T, 1, 0, MEM>(a_dx) {}

    inline double spectralRadius() const {return (2.0*DIM) / pow(this->dxMin(), 2); }
    inline void flux(
            BoxData<T, 1, MEM>& a_flux,
            const BoxData<T, 1, MEM>& a_state,
            int a_dir) const
    {
        auto GRAD = Stencil<T>::DiffCellToFace(a_dir, Side::Lo);
        //auto GRAD = 1.0*Shift::Zeros() - 1.0*Shift::Basis(a_dir, -1);
        a_flux |= GRAD(a_state, 1.0 / this->dx()[a_dir]); 
    }
};
