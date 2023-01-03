#pragma once
#include "Proto.H"

#define NUMCOMPS 1
#define NUMAUX 1

using namespace Proto;

template<typename T, MemType MEM = MEMTYPE_DEFAULT>
class BoxOp_TestSource : public BoxOp<T, NUMCOMPS, NUMAUX, MEM>
{
    public:

    // How many ghost cells does the operator need from the state variables
    inline static Point ghost() { return Point::Ones(1);}
    
    // How many ghost cells does the operator need from the auxiliary variables
    inline static Point auxGhost() { return Point::Ones(1);}
    
    // What is the intended order of accuracy of the operator
    inline static constexpr int order() { return 2; }

    inline BoxOp_TestSource() : BoxOp<T, NUMCOMPS, NUMAUX, MEM>() {}
    inline BoxOp_TestSource(double a_dx) : BoxOp<T, NUMCOMPS, NUMAUX, MEM>(a_dx) {}
    inline BoxOp_TestSource(Array<double, DIM> a_dx) : BoxOp<T, NUMCOMPS, NUMAUX, MEM>(a_dx)
    {
    }

    // The approximate spectral radius. Often needed if the operator is to be inverted using 
    // an iterative solver
    inline double spectralRadius() const {return 1;}
   
    inline void source(
            BoxData<T, NUMCOMPS>& a_source,
            const BoxData<T, NUMCOMPS>& a_state,
            const BoxData<T, NUMAUX>& a_aux) const
    {
        a_state.copyTo(a_source);
        a_source *= a_aux;
    }
};
