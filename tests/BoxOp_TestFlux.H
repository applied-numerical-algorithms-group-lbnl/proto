#pragma once
#include "Proto.H"

#define NUMCOMPS 1
#define NUMAUX 1

using namespace Proto;

template<typename T, MemType MEM = MEMTYPE_DEFAULT>
class BoxOp_TestFlux : public BoxOp<T, NUMCOMPS, NUMAUX, MEM>
{
    public:

    // How many ghost cells does the operator need from the state variables
    inline static Point ghost() { return Point::Ones(1);}
    
    // How many ghost cells does the operator need from the auxiliary variables
    inline static Point auxGhost() { return Point::Ones(1);}
    
    // What is the intended order of accuracy of the operator
    inline static constexpr int order() { return 2; }

    inline BoxOp_TestFlux() : BoxOp<T, NUMCOMPS, NUMAUX, MEM>() {}
    inline BoxOp_TestFlux(double a_dx) : BoxOp<T, NUMCOMPS, NUMAUX, MEM>(a_dx) {}
    inline BoxOp_TestFlux(Array<double, DIM> a_dx) : BoxOp<T, NUMCOMPS, NUMAUX, MEM>(a_dx)
    {
    }

    // The approximate spectral radius. Often needed if the operator is to be inverted using 
    // an iterative solver
    inline double spectralRadius() const {return 1;}
   
    inline void flux(
            BoxData<T, NUMCOMPS>& a_flux,
            const BoxData<T, NUMCOMPS>& a_state,
            const BoxData<T, NUMAUX>& a_aux,
            int a_dir) const
    {
        auto D = 1.0*Shift::Zeros() - 1.0*Shift::Basis(a_dir, -1);
        a_flux |= D(a_state, 1.0/this->dx()[a_dir]);
        a_flux *= a_aux;
    }
};
