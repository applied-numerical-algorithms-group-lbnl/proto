
#pragma once
#ifndef _PROTO_CROSS_PRODUCT_OP_
#define _PROTO_CROSS_PRODUCT_OP_

#include "Proto.H"

namespace Proto
{
    template <class T, template<class> typename PRODUCT_OP>
    class CrossProductOp : public CompositeOp<4>
    {

    protected:
        std::unique_ptr<PRODUCT_OP<T>> _productOp;
    public:
        inline CrossProductOp()
        {
        }

        inline std::array<Box, 4> domains(const Box &rangeBox) const
        {
            if (!_productOp)
            {
                std::array<Box, 4> domainBoxes;
                for (int ii = 0; ii < 4; ii++)
                {
                    domainBoxes[ii] = rangeBox;
                }
                return domainBoxes;
            }
            return _productOp->domains(rangeBox);
        }

        inline Box range(const std::array<Box, 4> &domainBoxes) const
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

        inline void
        operator()(
            BoxData<T,3>& AB_4thOrder,
            const BoxData<T,3>& A_4thOrder,
            const BoxData<T,3>& B_4thOrder,
            const BoxData<T,3>& A_2ndOrder,
            const BoxData<T,3>& B_2ndOrder)
        {
            Array<BoxData<T>,3> A2, B2, A4, B4;  
            for (int ii = 0; ii < 3; ii++)
            {
                A2[ii] = slice(A_2ndOrder, ii);
                A4[ii] = slice(A_4thOrder, ii);
                B2[ii] = slice(B_2ndOrder, ii);
                B4[ii] = slice(B_4thOrder, ii);
            }
            auto AB_01 = (*_productOp)(A4[0], B4[1], A2[0], B2[1]);
            auto AB_10 = (*_productOp)(A4[1], B4[0], A2[1], B2[0]);
            auto AB_12 = (*_productOp)(A4[1], B4[2], A2[1], B2[2]);
            auto AB_21 = (*_productOp)(A4[2], B4[1], A2[2], B2[1]);
            auto AB_20 = (*_productOp)(A4[2], B4[0], A2[2], B2[0]);
            auto AB_02 = (*_productOp)(A4[0], B4[2], A2[0], B2[2]);

            forallInPlace(
                []PROTO_LAMBDA(Var<T,3>& v_product,
                    Var<T>& v_01,Var<T>& v_10,
                    Var<T>& v_12,Var<T>& v_21,
                    Var<T>& v_20,Var<T>& v_02)
                {
                    v_product(0) = v_12(0) - v_21(0);
                    v_product(1) = v_20(0) - v_02(0);
                    v_product(2) = v_01(0) - v_10(0);
                }, AB_4thOrder, AB_01, AB_10, AB_12, AB_21, AB_20, AB_02);
        }

        inline BoxData<T, 3>
        operator()(
            const BoxData<T, 3> &A_4thOrder,
            const BoxData<T, 3> &B_4thOrder,
            const BoxData<T, 3> &A_2ndOrder,
            const BoxData<T, 3> &B_2ndOrder)
        {
            Box rangeBox = range(getDomains(A_4thOrder, B_4thOrder, A_2ndOrder, B_2ndOrder));
            BoxData<T, 3> output(rangeBox);
            (*this)(output, A_4thOrder, B_4thOrder, A_2ndOrder, B_2ndOrder);
            return output;
        }

        inline void
        operator()(
            LevelBoxData<T,3> &AB_4thOrder,
            const LevelBoxData<T,3> &A_4thOrder,
            const LevelBoxData<T,3> &B_4thOrder,
            const LevelBoxData<T,3> &A_2ndOrder,
            const LevelBoxData<T,3> &B_2ndOrder)
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