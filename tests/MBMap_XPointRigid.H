#pragma once

#include "Proto.H"

#define XPOINT_NUM_BLOCKS 4

using namespace Proto;

template<MemType MEM>
class MBMap_XPointRigid : public MBMapOp<MEM>
{
    public:

    using MBMapOp<MEM>::MBMapOp;

    //inline void setNumBlocks(unsigned int N) {m_numBlocks = N;}

    inline void apply(
            BoxData<double, DIM, MEM>& a_X,
            BoxData<double, 1, MEM>& a_J)
    {
        auto& dx = this->dx();
        auto X0 = this->X(a_X.box(), dx);
        auto block = this->block();
        int numBlocks = XPOINT_NUM_BLOCKS;
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
            }, a_X, X0, block, numBlocks
        );
        //a_J.setVal(sin(2.0*M_PI/numBlocks)*dx.product());
        a_J.setVal(sin(2.0*M_PI/numBlocks));
    }

    private:
};
