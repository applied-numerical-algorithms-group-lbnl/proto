
#pragma once
#ifndef _PROTO_PRODUCT_AVERAGE_OP_
#define _PROTO_PRODUCT_AVERAGE_OP_

#include "Proto.H"

namespace Proto
{
    template <class T>
    class ProductAverageOp : public ProductOrQuotientAverageOp<T, ProductAverageOp<T>>
    {

    private:

    public:

        template <unsigned int C, MemType MEM, unsigned int D = 1, unsigned int E = 1>
        inline void
        apply(
            BoxData<T, C, MEM, D, E> &AB_4thOrder,
            const BoxData<T, C, MEM, D, E> &A_4thOrder,
            const BoxData<T, C, MEM, D, E> &B_4thOrder,
            const BoxData<T, C, MEM, D, E> &A_2ndOrder,
            const BoxData<T, C, MEM, D, E> &B_2ndOrder)
        {
            PR_TIME("ProductAverageOp::apply");
            AB_4thOrder.setVal(0);
            for (auto& op : this->_diffOps)
            {
                BoxData<T,C,MEM,D,E> dA = op(A_2ndOrder);
                BoxData<T,C,MEM,D,E> dB = op(B_2ndOrder);
                AB_4thOrder.incrementProduct(dA, dB, 1.0/12.0);
            }
            AB_4thOrder.incrementProduct(A_4thOrder, B_4thOrder);
        }
    };
}
#endif