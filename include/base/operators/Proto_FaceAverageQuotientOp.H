
#pragma once
#ifndef _PROTO_FACE_AVERAGE_QUOTIENT_OP_
#define _PROTO_FACE_AVERAGE_QUOTIENT_OP_

#include "Proto.H"

namespace Proto
{
    template <class T>
    class FaceAverageQuotientOp : public QuotientAverageOp<T>
    {

    private:
        int _normDir;
        Array<Stencil<T>, DIM> _diffs;

    public:
        inline FaceAverageQuotientOp(int normDir)
        {
            for (int dd = 0; dd < DIM; dd++)
            {
                if (dd == normDir) { continue; }
                Stencil<T> diff = 0.5*Shift::Basis(dd, 1) - 0.5*Shift::Basis(dd, -1);
                this->_diffOps.push_back(diff);
            }
        }
    };
}
#endif