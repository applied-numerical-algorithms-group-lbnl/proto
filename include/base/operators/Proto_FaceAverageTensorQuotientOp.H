
#pragma once
#ifndef _PROTO_FACE_AVERAGE_TENSOR_QUOTIENT_OP_
#define _PROTO_FACE_AVERAGE_TENSOR_QUOTIENT_OP_

#include "Proto.H"

namespace Proto
{
    template <class T>
    class FaceAverageTensorQuotientOp : public CompositeOp<4>
    {

    private:
        int _normDir;
        std::unique_ptr<FaceAverageQuotientOp<T>> _quotient;

    public:
        inline FaceAverageTensorQuotientOp(int normDir) : _normDir(normDir)
        {
            _quotient = std::make_unique<FaceAverageQuotientOp<T>>(normDir);
        }

        inline std::array<Box, 4> domains(const Box &rangeBox) const
        {
            return _quotient->domains(rangeBox);
        }

        inline Box range(const std::array<Box, 4> &domainBoxes) const
        {
            return _quotient->range(domainBoxes);
        }

        template <unsigned int C, MemType MEM, unsigned int D = 1, unsigned int E = 1>
        inline void
        operator()(
            BoxData<T,C,MEM,D,E>& Q_4thOrder,
            const BoxData<T,C,MEM,D,E>& A_4thOrder,
            const BoxData<T,1,MEM,1,1>& B_4thOrder,
            const BoxData<T,C,MEM,D,E>& A_2ndOrder,
            const BoxData<T,1,MEM,1,1>& B_2ndOrder)  
        {
            PR_TIME("FaceAverageTensorQuotientOp::operator()");
            for (int ee = 0; ee < E;ee++)
            for (int dd = 0; dd < D; dd++)
            for (int cc = 0; cc < C; cc++)
            {
                auto Q = slice(Q_4thOrder,cc,dd,ee);
                auto A4 = slice(A_4thOrder,cc,dd,ee);
                auto A2 = slice(A_2ndOrder,cc,dd,ee);
                (*_quotient)(Q,A4,B_4thOrder,A2,B_2ndOrder);
            }
        }

        template <unsigned int C, MemType MEM, unsigned int D = 1, unsigned int E = 1>
        inline BoxData<T, C, MEM, D, E>
        operator()(
            const BoxData<T,C,MEM,D,E>& A_4thOrder,
            const BoxData<T,1,MEM,1,1>& B_4thOrder,
            const BoxData<T,C,MEM,D,E>& A_2ndOrder,
            const BoxData<T,1,MEM,1,1>& B_2ndOrder)  
        {
            Box rangeBox = range(getDomains(A_4thOrder, B_4thOrder, A_2ndOrder, B_2ndOrder));
            BoxData<T, C, MEM, D, E> output(rangeBox);
            (*this)(output, A_4thOrder, B_4thOrder, A_2ndOrder, B_2ndOrder);
            return output;
        }
    };
}
#endif