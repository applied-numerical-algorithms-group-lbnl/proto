#pragma once
#ifndef _BOXOP_ADVECTION_
#define _BOXOP_ADVECTION_
#include "ProtoAMR.H"

#define NUMCOMPS 1

using namespace Proto;

PROTO_KERNEL_START
void
f_thresholdF(const Point p, BoxData<short> &a_tags, BoxData<double,NUMCOMPS> &a_U)
{
  double thresh = .0001;
  if (a_U(p) > thresh) {a_tags(p) = 1;}
  else {a_tags(p) = 0;};
}
PROTO_KOKKOS_END(f_thresholdF, f_threshold);

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
     // forallInPlace(f_threshold, a_tags, a_state);
    }

    struct Advection_helper {
      Kokkos::View<short**> tags;
      Kokkos::View<T**> state;

      Advection_helper(short* tags_ptr, T* state_ptr, const Box box) : 
        tags(tags_ptr, box.size(0), box.size(1)), state(state_ptr, box.size(0), box.size(1)) {}

      KOKKOS_INLINE_FUNCTION
      void operator()(const int64_t i, const int64_t j) const {
          double thresh = .0001;
          if (state(i,j) > thresh) {tags(i,j) = 1;}
          else {tags(i,j) = 0;}
      }
    };

    static inline void kokkos_generateTags(
            TagData& a_tags, // BoxData<short,1>
            BoxData<T, NUMCOMPS>& a_state) // NUMCOMPS=1
    {
        Box box = a_state.box();
        std::tuple<TagData&,BoxData<T,NUMCOMPS>&> tuple(a_tags, a_state);
        Point lo(box.low()), hi(box.high());
        Kokkos::parallel_for("threshold", policy_t(
            {lo[0],lo[1]}, {hi[0],hi[1]}), struct_f_thresholdF<TagData&,BoxData<T,NUMCOMPS>&>(tuple));
    }
};
#endif //end include guard
