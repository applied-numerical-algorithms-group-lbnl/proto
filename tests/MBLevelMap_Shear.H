#pragma once

#include "Proto.H"

#define SHEAR_MAP_SLOPE 1.0

using namespace Proto;

template<MemType MEM>
class MBLevelMap_Shear : public MBLevelMap<MEM>
{
    public:

    using MBLevelMap<MEM>::MBLevelMap;

    inline void apply(
            BoxData<double, DIM, MEM>& a_X,
            BoxData<double, 1, MEM>& a_J,
            unsigned int a_block)
    {
        auto& layout = this->layout();
        auto& dx = this->dx(a_block);
        auto X0 = this->X(a_X.box(), dx);
        forallInPlace(
            [] PROTO_LAMBDA (Var<double, DIM, MEM>& a_X, Var<double, DIM, MEM>& a_X0, int a_block)
            {
                for (int dir = 2; dir < DIM; dir++)
                {
                    a_X(dir) = a_X0(dir);
                }
                double X = a_X0(0);
                double Y = a_X0(1);
                switch (a_block)
                {
                    // lower left
                    case 0:
                        a_X(0) = X - 1.0;
                        a_X(1) = Y - 1.0;
                        break;
                    // upper left
                    case 3:
                        a_X(0) = X - 1.0;
                        a_X(1) = Y;
                        break;
                    // lower right
                    case 1:
                        a_X(0) = X;
                        a_X(1) = Y + SHEAR_MAP_SLOPE*X - 1.0;
                        break;
                    // upper right
                    case 2:
                        a_X(0) = X;
                        a_X(1) = Y + SHEAR_MAP_SLOPE*X;
                        break;
                }
            }, a_X, X0, a_block
        );
        a_J.setVal(1.0);
    }
};

