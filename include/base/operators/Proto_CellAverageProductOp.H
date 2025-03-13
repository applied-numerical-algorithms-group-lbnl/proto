
#pragma once
#ifndef _PROTO_CELL_AVERAGE_PRODUCT_OP_
#define _PROTO_CELL_AVERAGE_PRODUCT_OP_

#include "Proto.H"

namespace Proto
{
    template <class T>
    class CellAverageProductOp : public CompositeOp
    {

    private:
        Array<Stencil<T>, DIM> _diffs;

    public:
        inline CellAverageProductOp()
        {
            for (int dd = 0; dd < DIM; dd++)
            {
                _diffs[dd] = 0.5*Shift::Basis(dd,1) - 0.5*Shift::Basis(dd,-1);
            }
        }

        inline std::vector<Box> domains(const Box &rangeBox) const
        {
            std::vector<Box> domainBoxes;

            domainBoxes.push_back(rangeBox);
            domainBoxes.push_back(rangeBox);

            Box diffBox;
            for (int dd = 0; dd < DIM; dd++)
            {
                diffBox += _diffs[dd].domain(rangeBox);
            }

            domainBoxes.push_back(diffBox);
            domainBoxes.push_back(diffBox);

            return domainBoxes;
        }

        inline Box range(const std::vector<Box> &domainBoxes) const
        {
            PROTO_ASSERT(domainBoxes.size() == 4,
                "CellAverageProduct::range | Error: Incorrect input size (%i). Should be 4.",
                domainBoxes.size());
            
            Box rangeBox = domainBoxes[0];
            rangeBox &= domainBoxes[1];
            for (int dd = 0; dd < DIM; dd++)
            {
                rangeBox &= _diffs[dd].range(domainBoxes[2]);
                rangeBox &= _diffs[dd].range(domainBoxes[3]);
            }
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
            PR_TIME("CellAverageProductOp::operator()");
            AB_4thOrder.setVal(0);
            for (int dd = 0; dd < DIM; dd++)
            {
                BoxData<T,C,MEM,D,E> dA = _diffs[dd](A_2ndOrder);
                BoxData<T,C,MEM,D,E> dB = _diffs[dd](B_2ndOrder);
                AB_4thOrder.incrementProduct(dA, dB, 1.0/12.0);
            }
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