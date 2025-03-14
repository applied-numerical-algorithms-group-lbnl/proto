
#pragma once
#ifndef _PROTO_MATRIX_PRODUCT_OP_
#define _PROTO_MATRIX_PRODUCT_OP_

#include "Proto.H"

namespace Proto
{
    template <class T, template<class> typename PRODUCT_OP>
    class MatrixProductOp : public CompositeOp
    {

    protected:
        std::shared_ptr<PRODUCT_OP<T>> _productOp;

    public:
        inline MatrixProductOp()
        {
        }

        inline std::vector<Box> domains(const Box &rangeBox) const
        {
            if (!_productOp)
            {
                std::vector<Box> domainBoxes;
                for (int ii = 0; ii < 4; ii++)
                {
                    domainBoxes.push_back(rangeBox);
                }
                return domainBoxes;
            }
            return _productOp->domains(rangeBox);
        }

        inline Box range(const std::vector<Box> &domainBoxes) const
        {
            if (!_productOp)
            {
                Box rangeBox = domainBoxes[0];
                for (auto box : domainBoxes)
                {
                    rangeBox &= box;
                }
                return rangeBox;
            }
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
                "FaceAverageMatrixProductOp::operator() | Error: Matrix inner dimensions must be equal (%u != %u)",
                DL, CR);
            PR_TIME("FaceAverageMatrixProductOp::operator()");
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
                if (!_productOp)
                {
                    AB.incrementProduct(A4, B4);
                } else {
                    AB += (*_productOp)(A4,B4,A2,B2);
                }
            }
        }

        template <unsigned int CL, unsigned int CR, MemType MEM,
            unsigned int DL, unsigned int DR, unsigned int E = 1>
        inline void
        transposeLeft(
            BoxData<T,DL,MEM,DR,E>& AB_4thOrder,
            const BoxData<T,CL,MEM,DL,E>& A_4thOrder,
            const BoxData<T,CR,MEM,DR,E>& B_4thOrder,
            const BoxData<T,CL,MEM,DL,E>& A_2ndOrder,
            const BoxData<T,CR,MEM,DR,E>& B_2ndOrder)
        {
            PROTO_ASSERT(CL == CR,
                "FaceAverageMatrixProductOp::transposeLeft | Error: Matrix inner dimensions must be equal (%u != %u)",
                CL, CR);
            PR_TIME("FaceAverageMatrixProductOp::transposeLeft");
            AB_4thOrder.setToZero();
            for (int ee = 0; ee < E;  ee++)
            for (int ii = 0; ii < DL; ii++)
            for (int jj = 0; jj < DR; jj++)
            for (int kk = 0; kk < CL; kk++)
            {
                auto AB = slice(AB_4thOrder, ii, jj, ee);
                const auto A4 = slice(A_4thOrder,  kk, ii, ee);
                const auto A2 = slice(A_2ndOrder,  kk, ii, ee);
                const auto B4 = slice(B_4thOrder,  kk, jj, ee);
                const auto B2 = slice(B_2ndOrder,  kk, jj, ee);
                if (!_productOp)
                {
                    AB.incrementProduct(A4, B4);
                } else {
                    AB += (*_productOp)(A4,B4,A2,B2);
                }
            }
        }

        template <unsigned int CL, unsigned int CR, MemType MEM,
            unsigned int DL, unsigned int DR, unsigned int E = 1>
        inline void
        transposeRight(
            BoxData<T,CL,MEM,CR,E>& AB_4thOrder,
            const BoxData<T,CL,MEM,DL,E>& A_4thOrder,
            const BoxData<T,CR,MEM,DR,E>& B_4thOrder,
            const BoxData<T,CL,MEM,DL,E>& A_2ndOrder,
            const BoxData<T,CR,MEM,DR,E>& B_2ndOrder)
        {
            PROTO_ASSERT(DL == DR,
                "FaceAverageMatrixProductOp::transposeLeft | Error: Matrix inner dimensions must be equal (%u != %u)",
                DL, DR);
            PR_TIME("FaceAverageMatrixProductOp::transposeLeft");
            AB_4thOrder.setToZero();
            for (int ee = 0; ee < E;  ee++)
            for (int ii = 0; ii < CL; ii++)
            for (int jj = 0; jj < CR; jj++)
            for (int kk = 0; kk < DL; kk++)
            {
                auto AB = slice(AB_4thOrder, ii, jj, ee);
                const auto A4 = slice(A_4thOrder,  ii, kk, ee);
                const auto A2 = slice(A_2ndOrder,  ii, kk, ee);
                const auto B4 = slice(B_4thOrder,  jj, kk, ee);
                const auto B2 = slice(B_2ndOrder,  jj, kk, ee);
                if (!_productOp)
                {
                    AB.incrementProduct(A4, B4);
                } else {
                    AB += (*_productOp)(A4,B4,A2,B2);
                }
            }
        }

        template <unsigned int CL, unsigned int CR, MemType MEM,
            unsigned int DL, unsigned int DR, unsigned int E = 1>
        inline BoxData<T, CL, MEM, DR, E>
        operator()(
            const BoxData<T,CL,MEM,DL,E>& A_4thOrder,
            const BoxData<T,CR,MEM,DR,E>& B_4thOrder,
            const BoxData<T,CL,MEM,DL,E>& A_2ndOrder,
            const BoxData<T,CR,MEM,DR,E>& B_2ndOrder)
        {
            Box rangeBox = range(getDomains(A_4thOrder, B_4thOrder, A_2ndOrder, B_2ndOrder));
            BoxData<T, CL, MEM, DR, E> output(rangeBox);
            (*this)(output, A_4thOrder, B_4thOrder, A_2ndOrder, B_2ndOrder);
            return output;
        }

        template <unsigned int CL, unsigned int CR, MemType MEM,
            unsigned int DL, unsigned int DR, unsigned int E = 1>
        inline BoxData<T, DL, MEM, DR, E>
        transposeLeft(
            const BoxData<T,CL,MEM,DL,E>& A_4thOrder,
            const BoxData<T,CR,MEM,DR,E>& B_4thOrder,
            const BoxData<T,CL,MEM,DL,E>& A_2ndOrder,
            const BoxData<T,CR,MEM,DR,E>& B_2ndOrder)
        {
            Box rangeBox = range(getDomains(A_4thOrder, B_4thOrder, A_2ndOrder, B_2ndOrder));
            BoxData<T, DL, MEM, DR, E> output(rangeBox);
            transposeLeft(output, A_4thOrder, B_4thOrder, A_2ndOrder, B_2ndOrder);
            return output;
        }

        template <unsigned int CL, unsigned int CR, MemType MEM,
            unsigned int DL, unsigned int DR, unsigned int E = 1>
        inline BoxData<T, CL, MEM, CR, E>
        transposeRight(
            const BoxData<T,CL,MEM,DL,E>& A_4thOrder,
            const BoxData<T,CR,MEM,DR,E>& B_4thOrder,
            const BoxData<T,CL,MEM,DL,E>& A_2ndOrder,
            const BoxData<T,CR,MEM,DR,E>& B_2ndOrder)
        {
            Box rangeBox = range(getDomains(A_4thOrder, B_4thOrder, A_2ndOrder, B_2ndOrder));
            BoxData<T, CL, MEM, CR, E> output(rangeBox);
            transposeRight(output, A_4thOrder, B_4thOrder, A_2ndOrder, B_2ndOrder);
            return output;
        }
    };
}
#endif