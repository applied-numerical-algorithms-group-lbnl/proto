
#pragma once
#ifndef _PROTO_EDGE_AVERAGE_PRODUCT_OP_
#define _PROTO_EDGE_AVERAGE_PRODUCT_OP_

#include "Proto.H"

namespace Proto
{
    template <class T>
    class EdgeAverageProductOp : public CompositeOp<4>
    {

    private:
        int _edgeDir;
        Stencil<T> _diff;

    public:
        inline EdgeAverageProductOp(int edgeDir) : _edgeDir(edgeDir)
        {
            _diff = 0.5*Shift::Basis(edgeDir, 1) - 0.5*Shift::Basis(edgeDir, -1);
        }

        inline std::array<Box, 4> domains(const Box &rangeBox) const
        {
            std::array<Box, 4> domainBoxes;
            Box diffBox = _diff.domain(rangeBox);
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
            rangeBox &= _diff.range(domainBoxes[2]);
            rangeBox &= _diff.range(domainBoxes[3]);
            return rangeBox;
        }

        template <unsigned int C, MemType MEM, unsigned int D = 1, unsigned int E = 1>
        inline void
        operator()(
            BoxData<T, C, MEM, D, E> &AB_4thOrder,
            const BoxData<T, C, MEM, D, E> &A_4thOrder,
            const BoxData<T, C, MEM, D, E> &B_4thOrder,
            const BoxData<T, C, MEM, D, E> &A_2ndOrder,
            const BoxData<T, C, MEM, D, E> &B_2ndOrder)
        {
            PR_TIME("EdgeAverageProductOp::operator()");
            AB_4thOrder.setVal(0);
            BoxData<T,C,MEM,D,E> dA = _diff(A_2ndOrder);
            BoxData<T,C,MEM,D,E> dB = _diff(B_2ndOrder);
            AB_4thOrder.incrementProduct(dA, dB, 1.0/12.0);
            AB_4thOrder.incrementProduct(A_4thOrder, B_4thOrder);
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