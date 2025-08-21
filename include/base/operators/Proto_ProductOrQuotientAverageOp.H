
#pragma once
#ifndef _PROTO_PRODUCT_OR_QUOTIENT_AVERAGE_OP_
#define _PROTO_PRODUCT_OR_QUOTIENT_AVERAGE_OP_

#include "Proto.H"

namespace Proto
{
    template <class T, typename OP>
    class ProductOrQuotientAverageOp : public CompositeOp<4>
    {

    protected:
        std::vector<Stencil<T>> _diffOps;

    public:
        inline ProductOrQuotientAverageOp()
        {
        }

        inline std::array<Box, 4> domains(const Box &rangeBox) const
        {
            std::array<Box, 4> domainBoxes;
            Box diffBox = rangeBox;
            for (auto& op : _diffOps)
            {
                diffBox += op.domain(rangeBox);
            }

            domainBoxes[0] = rangeBox;
            domainBoxes[1] = rangeBox;
            domainBoxes[2] = diffBox;
            domainBoxes[3] = diffBox;

            return domainBoxes;
        }

        inline Box range(const std::array<Box, 4> &domainBoxes) const
        {
            Box rangeBox = domainBoxes[0];
            rangeBox &= domainBoxes[1];
            for (auto& op : _diffOps)
            {
                rangeBox &= op.range(domainBoxes[2]);
                rangeBox &= op.range(domainBoxes[3]);
            }
            return rangeBox;
        }

        template <unsigned int C, MemType MEM, unsigned int D = 1, unsigned int E = 1>
        inline void
        operator()(
            BoxData<T,C,MEM,D,E>& AB_4thOrder,
            const BoxData<T, C, MEM, D, E> &A_4thOrder,
            const BoxData<T, C, MEM, D, E> &B_4thOrder,
            const BoxData<T, C, MEM, D, E> &A_2ndOrder,
            const BoxData<T, C, MEM, D, E> &B_2ndOrder)
        {
            static_cast<OP*>(this)->apply(AB_4thOrder, A_4thOrder, B_4thOrder, A_2ndOrder, B_2ndOrder);
        }

        template <unsigned int C, MemType MEM, unsigned int D = 1, unsigned int E = 1>
        inline BoxData<T, C, MEM, D, E>
        operator()(
            const BoxData<T, C, MEM, D, E> &A_4thOrder,
            const BoxData<T, C, MEM, D, E> &B_4thOrder,
            const BoxData<T, C, MEM, D, E> &A_2ndOrder,
            const BoxData<T, C, MEM, D, E> &B_2ndOrder)
        {
            Box rangeBox = range(getDomains(A_4thOrder, B_4thOrder, A_2ndOrder, B_2ndOrder));
            BoxData<T, C, MEM, D, E> output(rangeBox);
            (*this)(output, A_4thOrder, B_4thOrder, A_2ndOrder, B_2ndOrder);
            return output;
        }

        template <unsigned int C, MemType MEM, Centering CTR>
        inline void
        operator()(
            LevelBoxData<T, C, MEM, CTR> &AB_4thOrder,
            const LevelBoxData<T, C, MEM, CTR> &A_4thOrder,
            const LevelBoxData<T, C, MEM, CTR> &B_4thOrder,
            const LevelBoxData<T, C, MEM, CTR> &A_2ndOrder,
            const LevelBoxData<T, C, MEM, CTR> &B_2ndOrder)
        {
            PROTO_ASSERT(A_4thOrder.layout().compatible(B_4thOrder.layout()),
                "Incompatible layouts");
            PROTO_ASSERT(A_4thOrder.layout().compatible(A_2ndOrder.layout()),
                "Incompatible layouts");
            PROTO_ASSERT(A_4thOrder.layout().compatible(B_2ndOrder.layout()),
                "Incompatible layouts");
            PROTO_ASSERT(A_4thOrder.layout().compatible(AB_4thOrder.layout()),
                "Incompatible layouts");
            auto &layout = A_4thOrder.layout();
            for (auto index : layout)
            {
                const auto &A4 = A_4thOrder[index];
                const auto &B4 = B_4thOrder[index];
                const auto &A2 = A_2ndOrder[index];
                const auto &B2 = B_2ndOrder[index];
                auto AB4 = AB_4thOrder[index];
                (*this)(AB4, A4, B4, A2, B2);
            }
        }
    };
}
#endif