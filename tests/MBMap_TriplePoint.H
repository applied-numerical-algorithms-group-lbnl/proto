#pragma once

#include "Proto.H"

using namespace Proto;

namespace {
    ACCEL_DECORATION
    void f_triplePointMap(Array<double, DIM>& a_X, const Array<double, DIM>& a_X0, int a_block, int a_N)
    {
        double X = a_X0[0];
        double Y = a_X0[1];
#if DIM == 3
        double Z = a_X0[2];
#endif
        double dtheta_mapped = 2*M_PI/a_N;
        double dtheta_unmapped = M_PI/2.0;

        double u = (1-X);
        double v = (1-Y);
        double r = sqrt(u*u + v*v + 2*u*v*cos(dtheta_mapped));
        double alpha = dtheta_mapped;
        double theta = M_PI + atan2(v*sin(alpha), u + v*cos(alpha)) + alpha*a_block;

        a_X[0] = r*cos(theta);
        a_X[1] = r*sin(theta);
#if DIM == 3
        a_X[2] = Z;
#endif
    }
}

template<MemType MEM>
class MBMap_TriplePoint : public MBMapOp<MEM>
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
        int numBlocks = 3;
        forallInPlace(
            [] PROTO_LAMBDA (
                Var<double, DIM, MEM>& a_X, Var<double, DIM, MEM>& a_X0, int a_block, int a_N)
            {
                Array<double, DIM> X0;
                Array<double, DIM> X;
                for (int dir = 0; dir < DIM; dir++) { X0[dir] = a_X0(dir); };
                f_triplePointMap(X, X0, a_block, a_N);
                for (int dir = 0; dir < DIM; dir++) { a_X(dir) = X[dir]; }
            }, a_X, X0, block, numBlocks
        );
        a_J.setVal(sin(2.0*M_PI/numBlocks));
    }

    inline Array<double, DIM> apply(const Array<double, DIM>& a_X)
    {
        auto block = this->block();
        int numBlocks = 3;
        Array<double, DIM> X;
        f_triplePointMap(X,a_X, block, numBlocks);
        return X;
    }
    private:
};
