
#pragma once
#ifndef _PROTO_QUOTIENT_AVERAGE_OP_
#define _PROTO_QUOTIENT_AVERAGE_OP_

#include "Proto.H"

namespace Proto
{
    template <class T>
    class QuotientAverageOp : public ProductOrQuotientAverageOp<T, QuotientAverageOp<T>>
    {

    public:

    template <unsigned int C, MemType MEM, unsigned int D = 1, unsigned int E = 1>
    inline void
    apply(
        BoxData<T, C, MEM, D, E> &Quotient_4thOrder,
        const BoxData<T, C, MEM, D, E> &A_4thOrder,
        const BoxData<T, C, MEM, D, E> &B_4thOrder,
        const BoxData<T, C, MEM, D, E> &A_2ndOrder,
        const BoxData<T, C, MEM, D, E> &B_2ndOrder)
    {
        auto& Q = Quotient_4thOrder;
        PR_TIME("QuotientAverageOp::operator()");
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
        for (auto& op : this->_diffOps)
        {
            BoxData<T,C,MEM,D,E> dA = op(A_2ndOrder);
            BoxData<T,C,MEM,D,E> dB = op(B_2ndOrder);
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
    };
}
#endif