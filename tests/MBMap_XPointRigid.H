#pragma once

#include "Proto.H"

using namespace Proto;

class MBMap_XPointRigid : public MBMapOp
{
    public:

    using MBMapOp::MBMapOp;

    inline void setNumBlocks(unsigned int N) {m_numBlocks = N;}

    template<MemType MEM>
    inline void apply(
            BoxData<double, DIM, MEM>& a_X,
            BoxData<double, 1, MEM>& a_J,
            unsigned int a_block)
    {
        auto& layout = this->layout();
        auto& dx = this->dx(a_block);
        auto X0 = this->X(a_X.box(), dx);
        forallInPlace(
            [] PROTO_LAMBDA (
                Var<double, DIM, MEM>& a_X, Var<double, DIM, MEM>& a_X0, int a_block, int a_N)
            {
                double X = a_X0(0);
                double Y = a_X0(1);
#if DIM == 3
                double Z = a_X0(2);
#endif
                double dtheta_mapped = 2*M_PI/a_N;
                double dtheta_unmapped = M_PI/2.0;
                double theta_ratio = dtheta_mapped/dtheta_unmapped;

                double u = (1-X);
                double v = (1-Y);
                double r = sqrt(u*u + v*v + 2*u*v*cos(dtheta_mapped));
                double alpha = dtheta_mapped;
                double theta = M_PI + atan2(v*sin(alpha), u + v*cos(alpha)) + alpha*a_block;

                a_X(0) = r*cos(theta);
                a_X(1) = r*sin(theta);
#if DIM == 3
                a_X(2) = Z;
#endif
            }, a_X, X0, a_block, m_numBlocks
        );
        a_J.setVal(sin(2.0*M_PI/m_numBlocks));
    }

    private:

    unsigned int m_numBlocks;
};
