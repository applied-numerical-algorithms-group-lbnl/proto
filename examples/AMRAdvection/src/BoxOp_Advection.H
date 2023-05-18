#pragma once
#ifndef _BOXOP_ADVECTION_
#define _BOXOP_ADVECTION_
#include "ProtoAMR.H"

#define NUMCOMPS 1

using namespace Proto;

#ifdef PROTO_ACCEL
PROTO_KERNEL_START
void
k_threshold(const Point p, DevData<short> &a_tags, DevData<double,NUMCOMPS> &a_U)
{
  double thresh = .0001;
  if (a_U(p) > thresh) {a_tags(p) = 1;}
  else {a_tags(p) = 0;};
}
PROTO_KOKKOS_END(k_thresholdF, k_threshold);
#endif
PROTO_KERNEL_START
void
f_thresholdF(Var<short> &a_tags, const Var<double,NUMCOMPS> &a_U)
{
  double thresh = .0001;
  if (a_U(0) > thresh) {a_tags(0) = 1;}
  else {a_tags(0) = 0;};
}
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
            TagData& a_tags, // BoxData<short,1>
            BoxData<T, NUMCOMPS>& a_state) {
      forallInPlace(f_threshold, a_tags, a_state);
    }

#ifdef PROTO_ACCEL
    static inline void kokkos_generateTags(
            TagData& a_tags, // BoxData<short,1>
            BoxData<T, NUMCOMPS>& a_state) // NUMCOMPS=1
    {
        DevData<short> d_tags(a_tags); DevData<double, NUMCOMPS> d_state(a_state);
        Kokkos::parallel_for("threshold", policy(a_tags.box()&a_state.box()), k_thresholdF(d_tags,d_state));
    }
#endif
};
#endif //end include guard
