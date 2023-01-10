#pragma once
#ifndef _PROTO_INTERP_STENCIL_
#define _PROTO_INTERP_STENCIL_

#include "Proto_InterpStencil1D.H"

namespace Proto
{
    template<typename T>
    class InterpStencil
    {
        public:

        inline static InterpStencil<T> Constant(Point a_refRatio);
        inline static InterpStencil<T> Linear(Point a_refRatio);
        inline static InterpStencil<T> Quadratic(Point a_refRatio);
        inline static InterpStencil<T> Cubic(Point a_refRatio);
        inline static InterpStencil<T> Quartic(Point a_refRatio);

        InterpStencil(){};
        InterpStencil(Point a_refRatio);
        void define(Point a_refRatio);
        
        template<unsigned int C, MemType MEM=MEMTYPE_DEFAULT, unsigned char D, unsigned char E>
        inline void apply(BoxData<T,C,MEM,D,E>& a_output, BoxData<T,C,MEM,D,E>& a_input, T a_scale = 1.0);
        inline Point ghost() const;
        inline Box span() const;

        private:

        Stencil<T> build1D5pt(Array<double,5>& a_coefs, int a_dir);

        Array<InterpStencil1D<T>, DIM> m_interp;

    }; //end class InterpStencil
#include "implem/Proto_InterpStencilImplem.H"
#include "implem/Proto_InterpStencilDefs.H"
} // end Proto namespace
#endif //end include guard
