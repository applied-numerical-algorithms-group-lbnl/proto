
#pragma once
#ifndef _PROTO_FACE_AVERAGE_QUOTIENT_OP_
#define _PROTO_FACE_AVERAGE_QUOTIENT_OP_

#include "Proto.H"

namespace Proto
{
    template <class T>
    class FaceAverageQuotientOp : public CompositeOp
    {

    private:
        int _normDir;
        Array<Stencil<T>, DIM> _diffs;

    public:
        inline FaceAverageQuotientOp(int normDir)
        {
            _normDir = normDir;
            for (int dd = 0; dd < DIM; dd++)
            {
                if (dd == normDir)
                {
                    _diffs[dd] = 1.0*Shift::Zeros();
                } else {
                    _diffs[dd] = 0.5*Shift::Basis(dd, 1) - 0.5*Shift::Basis(dd, -1);
                }
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
                "FaceAverageQuotient::range | Error: Incorrect input size (%i). Should be 4.",
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
            BoxData<T, C, MEM, D, E> &Quotient_4thOrder,
            const BoxData<T, C, MEM, D, E> &A_4thOrder,
            const BoxData<T, C, MEM, D, E> &B_4thOrder,
            const BoxData<T, C, MEM, D, E> &A_2ndOrder,
            const BoxData<T, C, MEM, D, E> &B_2ndOrder)
        {
            auto& Q = Quotient_4thOrder;
            PR_TIME("CellAverageQuotientOp::operator()");
            forallInPlace(
                [=] PROTO_LAMBDA (
                    Var<T,C,MEM,D,E>& v_Q,
                    Var<T,C,MEM,D,E>& v_A,
                    Var<T,C,MEM,D,E>& v_B)
                {
                    for (int cc = 0; cc < C; cc++)
                    for (int dd = 0; dd < D; dd++)
                    for (int ee = 0; ee < E; ee++)
                    {
                        v_Q(cc,dd,ee) = (v_A(cc,dd,ee) / v_B(cc,dd,ee));
                    }
                }, Q, A_4thOrder, B_4thOrder);
            
            BoxData<T,C,MEM,D,E> dBSquaredFactor(Q.box());
            BoxData<T,C,MEM,D,E> dAdBTerm(Q.box());
            
            dBSquaredFactor.setVal(1.0);
            dAdBTerm.setToZero();
            for (int dd = 0; dd < DIM; dd++)
            {
                if (dd == _normDir) { continue; }
                BoxData<T,C,MEM,D,E> dA = _diffs[dd](A_2ndOrder);
                BoxData<T,C,MEM,D,E> dB = _diffs[dd](B_2ndOrder);
                forallInPlace(
                    [] PROTO_LAMBDA (
                        Var<T,C,MEM,D,E>& v_dBSquared,
                        Var<T,C,MEM,D,E>& v_dAdBTerm,
                        Var<T,C,MEM,D,E>& v_dA,
                        Var<T,C,MEM,D,E>& v_dB,
                        Var<T,C,MEM,D,E>& v_B)
                    {
                        for (int cc = 0; cc < C; cc++)
                        for (int dd = 0; dd < D; dd++)
                        for (int ee = 0; ee < E; ee++)
                        {
                            T BB = v_B(cc,dd,ee)*v_B(cc,dd,ee)*12.0;
                            v_dBSquared(cc,dd,ee) += v_dB(cc,dd,ee)*v_dB(cc,dd,ee)/BB;
                            v_dAdBTerm(cc,dd,ee) -= v_dA(cc,dd,ee)*v_dB(cc,dd,ee)/BB;
                        }
                    }, dBSquaredFactor, dAdBTerm, dA, dB, B_4thOrder);
            }
            Q *= dBSquaredFactor;
            Q += dAdBTerm;
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
            LevelBoxData<T, C, MEM, CTR> &Quotient_4thOrder,
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
            PROTO_ASSERT(A_4thOrder.layout().compatible(Quotient_4thOrder.layout()),
                "Incompatible layouts");
            auto &layout = A_4thOrder.layout();
            for (auto index : layout)
            {
                const auto &A4 = A_4thOrder[index];
                const auto &B4 = B_4thOrder[index];
                const auto &A2 = A_2ndOrder[index];
                const auto &B2 = B_2ndOrder[index];
                auto Q4 = Quotient_4thOrder[index];
                (*this)(Q4, A4, B4, A2, B2);
            }
        }
    };
}
#endif