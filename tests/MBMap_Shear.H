#pragma once

#include "Proto.H"

#define SHEAR_MAP_SLOPE 1.0

using namespace Proto;

namespace {
    ACCEL_DECORATION
    void f_shearMap(Array<double, DIM>& a_X, const Array<double, DIM>& a_X0, int a_block)
    {
        for (int dir = 2; dir < DIM; dir++)
        {
            a_X[dir] = a_X0[dir];
        }
        double X = a_X0[0];
        double Y = a_X0[1];
        switch (a_block)
        {
            // lower left
            case 0:
                a_X[0] = X - 1.0;
                a_X[1] = Y - 1.0;
                break;
                // upper left
            case 3:
                a_X[0] = X - 1.0;
                a_X[1] = Y;
                break;
                // lower right
            case 1:
                a_X[0] = X;
                a_X[1] = Y + SHEAR_MAP_SLOPE*X - 1.0;
                break;
                // upper right
            case 2:
                a_X[0] = X;
                a_X[1] = Y + SHEAR_MAP_SLOPE*X;
                break;
        }
    }
}

template<MemType MEM>
class MBMap_Shear : public MBMapOp<MEM>
{
    public:

    using MBMapOp<MEM>::MBMapOp;

    inline void apply(
            BoxData<double, DIM, MEM>& a_X,
            BoxData<double, 1, MEM>& a_J)
    {
        auto& dx = this->dx();
        auto X0 = this->X(a_X.box(), dx);
        auto block = this->block();
        forallInPlace(
        [] PROTO_LAMBDA (Var<double, DIM, MEM>& a_X, Var<double, DIM, MEM>& a_X0, int a_block)
        {
            Array<double, DIM> X0;
            for (int ii = 0; ii < DIM; ii++) { X0[ii] = a_X0(ii); }
            Array<double, DIM> X;
            f_shearMap(X, X0, a_block);
            for (int ii = 0; ii < DIM; ii++) { a_X(ii) = X[ii]; }
        }, a_X, X0, block
        );
        a_J.setVal(1.0);
    }
    
    inline Array<double, DIM> apply(const Array<double, DIM>& a_X)
    {
        auto block = this->block();
        Array<double, DIM> X;
        f_shearMap(X, a_X, block);
        return X;
    }
};

