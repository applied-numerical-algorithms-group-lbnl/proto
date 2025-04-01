
#pragma once
#ifndef _PROTO_EDGE_AVERAGE_PRODUCT_OP_
#define _PROTO_EDGE_AVERAGE_PRODUCT_OP_

#include "Proto.H"

namespace Proto
{
    template <class T>
    class EdgeAverageProductOp : public ProductAverageOp<T>
    {
    public:
        inline EdgeAverageProductOp(int edgeDir)
        {
            Stencil<T> diff = 0.5*Shift::Basis(edgeDir, 1) - 0.5*Shift::Basis(edgeDir, -1);
            this->_diffOps.push_back(diff);
        }
    };
}
#endif