#pragma once
#ifndef _PROTO_INTERP_STENCIL_
#define _PROTO_INTERP_STENCIL_

#include "Proto_InterpStencil1D.H"

namespace Proto
{
    template<typename T, unsigned int C, MemType MEM,
        unsigned char D, unsigned char E>
    struct LazyInterpStencil;

    template<typename T>
    class InterpStencil
    {
        public:

        inline static InterpStencil<T> Constant(Point a_refRatio);
        inline static InterpStencil<T> Linear(Point a_refRatio);
        inline static InterpStencil<T> Quadratic(Point a_refRatio);
        inline static InterpStencil<T> FiniteVolume(Point a_refRatio, int a_order = 5);

        InterpStencil(){};
        InterpStencil(Point a_refRatio);
        void define(Point a_refRatio);
        
        inline Box domain(Box a_range) const;
        inline Box range(Box a_domain) const;
        inline Point ghost() const;
        inline Box span() const;
        inline Point ratio() const {return m_ratio; } 

        template<unsigned int C, MemType MEM=MEMTYPE_DEFAULT, unsigned char D, unsigned char E>
        inline void apply(
                BoxData<T,C,MEM,D,E>& a_output,
                const BoxData<T,C,MEM,D,E>& a_input,
                Box a_box,
                bool a_overwrite,
                T a_scale = 1.0) const;

        template<unsigned int C, MemType MEM, unsigned char D, unsigned char E>
        inline LazyInterpStencil<T, C, MEM, D, E> operator()(
                const BoxData<T, C, MEM, D, E>& a_input,
                T a_scale = 1.0) const;

        template<unsigned int C, MemType MEM, unsigned char D, unsigned char E>
        inline LazyInterpStencil<T, C, MEM, D, E> operator()(
                const BoxData<T, C, MEM, D, E>& a_input,
                Box a_range,
                T a_scale = 1.0) const;

        private:
        
        Point m_ratio;
        Array<InterpStencil1D<T>, DIM> m_interp;

    }; //end class InterpStencil
    
    template<typename T, unsigned int C, MemType MEM,
        unsigned char D, unsigned char E>
    struct LazyInterpStencil
    {
        inline LazyInterpStencil(){}
        inline LazyInterpStencil(
                const InterpStencil<T>& a_stencil,
                const BoxData<T, C, MEM, D, E>& a_src,
                Box a_range, T a_scale);
        
        inline void apply(BoxData<T,C,MEM,D,E>& a_dst, bool a_overwrite) const;
        inline Box inferredRange() const;
        Box range;
        const InterpStencil<T>* op;
        const BoxData<T, C, MEM, D, E>* srcData;
        T scale;
    }; // end class LazyInterpStencil

#include "implem/Proto_InterpStencilImplem.H"
#include "implem/Proto_InterpStencilDefs.H"
   
    
    // Increment operator
    template<typename T, unsigned int C, MemType MEM,
        unsigned char D, unsigned char E>
    BoxData<T,C,MEM,D,E>& operator += (
            BoxData<T,C,MEM,D,E>& a_dst,
            const LazyInterpStencil<T,C,MEM,D,E>&& a_op)
    {
        a_op.apply(a_dst, false);
        return a_dst;
    }
    // Replace operator 
    template<typename T, unsigned int C, MemType MEM,
        unsigned char D, unsigned char E>
    BoxData<T,C,MEM,D,E>& operator |= (
            BoxData<T,C,MEM,D,E>& a_dst,
            const LazyInterpStencil<T,C,MEM,D,E>&& a_op)
    {
        a_op.apply(a_dst, true);
        return a_dst;
    }
} // end Proto namespace
#endif //end include guard
