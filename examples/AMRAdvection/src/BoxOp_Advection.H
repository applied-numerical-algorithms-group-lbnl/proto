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
  double thresh = .0001;
  if (a_U(0) > thresh) {a_tags(0) = 1;}
  else {a_tags(0) = 0;};
};

PROTO_KERNEL_END(f_thresholdF, f_threshold);

template<typename T, MemType MEM = MEMTYPE_DEFAULT>
class BoxOp_Advection : public BoxOp<T, NUMCOMPS, 0, MEM>
{
    public:

    using BoxOp<T,NUMCOMPS,0,MEM>::BoxOp;

    // These functions are not optional
    inline static Point ghost() { return Point::Ones(3);}
    inline static Point auxGhost() { return Point::Zeros();}
    inline static constexpr int order() { return 4; }

    //inline T spectralRadius() const {return (2.0*DIM) / pow(this->dxMin(), 2); }
  
    // Apply BCs by filling ghost cells in stage values. For Euler, this is done by calling
    // exchange. For the MHD code, it will be more complicated.
    // The interface is very provisional. We expect it to evolve as we d more real problems.
    inline void bcStage(
                        LevelBoxData<T,NUMCOMPS>& a_UStage,
                        const LevelBoxData<T,NUMCOMPS>& a_U0,
                        int a_stage)
    {
      a_UStage.exchange();
    }
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
