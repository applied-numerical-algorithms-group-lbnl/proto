#pragma once

#include "Proto.H"
//#define XPOINT_NBLOCKS 4
using namespace Proto;

namespace {
    ACCEL_DECORATION
    void f_xpointMap(Array<double, DIM>& a_X, const Array<double, DIM>& a_X0, int a_block, int a_N)
    {
        double X = a_X0[0];
        double Y = a_X0[1];

        double dtheta_mapped = 2*M_PI/a_N;
        double dtheta_unmapped = M_PI/2.0;

        double u = (1-X);
        double v = (1-Y);
        double r = sqrt(u*u + v*v + 2*u*v*cos(dtheta_mapped));
        double alpha = dtheta_mapped;
        double theta = M_PI + atan2(v*sin(alpha), u + v*cos(alpha)) + alpha*a_block;

        a_X[0] = r*cos(theta);
        a_X[1] = r*sin(theta);

        for (int ii = 2; ii < DIM; ii++)
        {
            a_X[ii] = a_X0[ii];
        }
    }

    ACCEL_DECORATION
    void f_inverseXPointMap(Array<double, DIM>& a_X0, const Array<double, DIM>& a_X, int a_block, int a_N)
    {

        double x = a_X[0];
        double y = a_X[1];
        
        double r = sqrt(x*x + y*y);
        double theta = atan2(y, x);
        
        double alpha = 2*M_PI/a_N;
        
        double gamma = theta - M_PI - alpha*a_block;
        double v = r*sin(gamma)/sin(alpha);
        double u = r*cos(gamma) - v*cos(alpha);
        
        a_X0[0] = 1.0 - u;
        a_X0[1] = 1.0 - v;
        
        for (int ii = 2; ii < DIM; ii++)
        {
            a_X0[ii] = a_X[ii];
        }

    }
}

template<int NBLOCK, MemType MEM>
class MBMap_XPointRigid : public MBMapOp<MEM>
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
        int n = NBLOCK;
        forallInPlace(
            [] PROTO_LAMBDA (
                Var<double, DIM, MEM>& a_X, Var<double, DIM, MEM>& a_X0, int a_block, int a_N)
            {
                Array<double, DIM> X0;
                Array<double, DIM> X;
                for (int dir = 0; dir < DIM; dir++) { X0[dir] = a_X0(dir); };
                f_xpointMap(X, X0, a_block, a_N);
                for (int dir = 0; dir < DIM; dir++) { a_X(dir) = X[dir]; }
            }, a_X, X0, block, n
        );
        a_J.setVal(sin(2.0*M_PI/n));
    }

    inline Array<double, DIM> apply(const Array<double, DIM>& a_X)
    {
        auto block = this->block();
        Array<double, DIM> X;
        f_xpointMap(X,a_X, block, NBLOCK);
        return X;
    }

    inline void inverse(
        BoxData<double, DIM, MEM>& a_X0,
        const BoxData<double, DIM, MEM>& a_X)
    {
        auto block = this->block();
        forallInPlace(
        [] PROTO_LAMBDA (Var<double, DIM, MEM>& a_X0, Var<double, DIM, MEM>& a_X, int a_block)
        {
            Array<double, DIM> X;
            for (int ii = 0; ii < DIM; ii++) { X[ii] = a_X(ii); }
            Array<double, DIM> X0;
            f_inverseXPointMap(X0, X, a_block, NBLOCK);
            for (int ii = 0; ii < DIM; ii++) { a_X0(ii) = X0[ii]; }
        }, a_X0, a_X, block
        );
    }
    
    inline Array<double, DIM> inverse(const Array<double, DIM>& a_X)
    {
        auto block = this->block();
        Array<double, DIM> X0;
        f_inverseXPointMap(X0, a_X, block, NBLOCK);
        return X0;
    }
    private:
};
