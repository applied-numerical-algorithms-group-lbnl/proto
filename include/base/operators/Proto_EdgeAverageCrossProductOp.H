
#pragma once
#ifndef _PROTO_EDGE_AVERAGE_CROSS_PRODUCT_OP_
#define _PROTO_EDGE_AVERAGE_CROSS_PRODUCT_OP_

#include "Proto.H"

namespace Proto
{
    template <class T>
    class EdgeAverageCrossProductOp : public CrossProductOp<T, EdgeAverageProductOp>
    {

    public:
        inline EdgeAverageCrossProductOp(int edgeDir)
        {
            this->_productOp = std::make_unique<EdgeAverageProductOp<T>>(edgeDir);
        }
    };
}
#endif