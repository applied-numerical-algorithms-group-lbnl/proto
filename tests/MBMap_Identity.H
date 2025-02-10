#pragma once

#include "Proto.H"

using namespace Proto;

template<MemType MEM>
class MBMap_Identity : public MBMapOp<MEM>
{
    public:

    using MBMapOp<MEM>::MBMapOp;

    inline void apply(
            BoxData<double, DIM, MEM>& a_X,
            BoxData<double, 1, MEM>& a_J)
    {
        auto& dx = this->dx();
        auto X0 = this->X(a_X.box(), dx);
        X0.copyTo(a_X);
        a_J.setVal(1.0);
    }

    inline Array<double, DIM> apply(const Array<double, DIM>& a_X)
    {
        return a_X;
    }
    private:
};
