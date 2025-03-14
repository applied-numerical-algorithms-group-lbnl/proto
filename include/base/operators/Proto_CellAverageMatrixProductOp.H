
#pragma once
#ifndef _PROTO_CELL_AVERAGE_MATRIX_PRODUCT_OP_
#define _PROTO_CELL_AVERAGE_MATRIX_PRODUCT_OP_

#include "Proto.H"

namespace Proto
{
    template <class T>
    class CellAverageMatrixProductOp : public CompositeOp
    {

    private:
        std::shared_ptr<CellAverageProductOp<T>> _productOp;

    public:
        inline CellAverageMatrixProductOp()
        {
            _productOp = std::make_shared<CellAverageProductOp<T>>();
        }

        inline std::vector<Box> domains(const Box &rangeBox) const
        {
            return _productOp->domains(rangeBox);
        }

        inline Box range(const std::vector<Box> &domainBoxes) const
        {
            return _productOp->range(domainBoxes);
        }

        template <unsigned int CL, unsigned int CR, MemType MEM,
            unsigned int DL, unsigned int DR, unsigned int E = 1>
        inline void
        operator()(
            BoxData<T,CL,MEM,DR,E>& AB_4thOrder,
            const BoxData<T,CL,MEM,DL,E>& A_4thOrder,
            const BoxData<T,CR,MEM,DR,E>& B_4thOrder,
            const BoxData<T,CL,MEM,DL,E>& A_2ndOrder,
            const BoxData<T,CR,MEM,DR,E>& B_2ndOrder)
        {
            PROTO_ASSERT(DL == CR,
                "CellAverageMatrixProductOp::operator() | Error: Matrix inner dimensions must be equal (%u != %u)",
                CR, DL);
            PR_TIME("CellAverageMatrixProductOp::operator()");
            AB_4thOrder.setToZero();
            for (int ee = 0; ee < E;  ee++)
            for (int ii = 0; ii < CL; ii++)
            for (int jj = 0; jj < DR; jj++)
            for (int kk = 0; kk < DL; kk++)
            {
                auto AB = slice(AB_4thOrder, ii, jj, ee);
                const auto B4 = slice(B_4thOrder,  kk, jj, ee);
                const auto B2 = slice(B_2ndOrder,  kk, jj, ee);
                const auto A4 = slice(A_4thOrder,  ii, kk, ee);
                const auto A2 = slice(A_2ndOrder,  ii, kk, ee);
                AB += (*_productOp)(A4,B4,A2,B2);
            }
        }

        // Optional, but easy to implement
        // template <unsigned int C, MemType MEM, unsigned int D = 1, unsigned int E = 1>
        // inline BoxData<T, C, MEM, D, E>
        // operator()(
        //     const BoxData<T, C, MEM, D, E> &input)
        // {
        //     Box rangeBox = range(getDomains(input));
        //     BoxData<T, C, MEM, D, E> output(rangeBox);
        //     (*this)(output, input);
        //     return output;
        // }
    };
}
#endif