#pragma once
#ifndef _BOXOP_ADVECTION_
#define _BOXOP_ADVECTION_
#include "ProtoAMR.H"

#define NUMCOMPS 1

using namespace Proto;

PROTO_KERNEL_START
void
f_thresholdF(
           Var<short>& a_tags,
           Var<double, NUMCOMPS>& a_U)
{
    double thresh = 1.0001;
    if (a_U(0) > thresh) {a_tags(0) = 1;}
    else {a_tags(0) = 0;};
};

PROTO_KERNEL_END(f_thresholdF, f_threshold);

template<typename T, MemType MEM = MEMTYPE_DEFAULT>
class BoxOp_Advection : public BoxOp<T, NUMCOMPS, 0, MEM>
{
    public:

    // These functions are not optional
    inline static Point ghost() { return Point::Ones(3);}
    inline static Point auxGhost() { return Point::Zeros();}
    inline static constexpr int order() { return 4; }

    inline BoxOp_Advection() : BoxOp<T, 1, 0, MEM>() {}
    inline BoxOp_Advection(T a_dx) : BoxOp<T, 1, 0, MEM>(a_dx) {}
    inline BoxOp_Advection(Array<T, DIM> a_dx) : BoxOp<T, 1, 0, MEM>(a_dx) {}

    //inline T spectralRadius() const {return (2.0*DIM) / pow(this->dxMin(), 2); }
    inline void flux(
            BoxData<T, NUMCOMPS>& a_flux,
            const BoxData<T, NUMCOMPS>& a_state,
            int a_dir) const
    {
        T vel = 1.0;
        Stencil<T> interp;
        if (vel > 0)
        {
            interp = Stencil<double>::CellToFace(a_dir,Side::Lo,5);
        } else {
            interp = Stencil<double>::CellToFace(a_dir,Side::Hi,5)*((1.0)*Shift::Basis(a_dir,-1));
        }
        a_flux |= interp(a_state, -vel);
    }
    
    static inline void generateTags(
            TagData& a_tags,
            BoxData<T, NUMCOMPS>& a_state)
    {
        forallInPlace(f_threshold, a_tags, a_state);
    }
};
#endif //end include guard
