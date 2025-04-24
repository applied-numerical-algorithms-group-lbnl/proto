
#pragma once
#ifndef _PROTO_CONVOLVE_FACE_OP_
#define _PROTO_CONVOLVE_FACE_OP_

namespace Proto
{

    template <class T>
    class ConvolveFaceOp : public CompositeOp<2>
    {

    private:
        Stencil<T> m_Lf;

    public:
        inline ConvolveFaceOp(int dir)
        {
            m_Lf = Stencil<T>::LaplacianFace(dir);
        }

        inline std::array<Box, 2> domains(const Box &rangeBox) const
        {
            std::array<Box, 2> domainBoxes;
            domainBoxes[0] = rangeBox;
            domainBoxes[1] = m_Lf.domain(rangeBox);
            return domainBoxes;
        }

        inline Box range(const std::array<Box, 2> &domainBoxes) const
        {
            return (m_Lf.range(domainBoxes[1]) & domainBoxes[0]);
        }

        template <unsigned int C, MemType MEM, unsigned int D = 1, unsigned int E = 1>
        inline void
        operator()(
            BoxData<T, C, MEM, D, E> &cellAveraged_4thOrder,
            const BoxData<T, C, MEM, D, E> &cellCentered_4thOrder,
            const BoxData<T, C, MEM, D, E> &cellCentered_2ndOrder)
        {
            PR_TIME("ConvolveFaceOp::operator()");
            cellAveraged_4thOrder |= m_Lf(cellCentered_2ndOrder, 1.0 / 24.0);
            cellAveraged_4thOrder += cellCentered_4thOrder;
        }

        template <unsigned int C, MemType MEM, unsigned int D = 1, unsigned int E = 1>
        inline BoxData<T, C, MEM, D, E>
        operator()(
            const BoxData<T, C, MEM, D, E> &cellCentered_4thOrder,
            const BoxData<T, C, MEM, D, E> &cellCentered_2ndOrder)
        {
            Box rangeBox = range(getDomains(cellCentered_4thOrder, cellCentered_2ndOrder));
            BoxData<T, C, MEM, D, E> cellAveraged(rangeBox);
            (*this)(cellAveraged, cellCentered_4thOrder, cellCentered_2ndOrder);
            return cellAveraged;
        }

        template <unsigned int C, MemType MEM, Centering CTR>
        inline void
        operator()(
            LevelBoxData<T, C, MEM, CTR> &cellAveraged_4thOrder,
            const LevelBoxData<T, C, MEM, CTR> &cellCentered_4thOrder,
            const LevelBoxData<T, C, MEM, CTR> &cellCentered_2ndOrder)
        {
            PROTO_ASSERT(cellAveraged_4thOrder.layout().compatible(cellCentered_4thOrder.layout()),
                         "incompatible layouts");
            PROTO_ASSERT(cellAveraged_4thOrder.layout().compatible(cellCentered_2ndOrder.layout()),
                         "incompatible layouts");

            auto &layout = cellAveraged_4thOrder.layout();
            for (auto index : layout)
            {
                auto &avg = cellAveraged_4thOrder[index];
                const auto &ctr4 = cellCentered_4thOrder[index];
                const auto &ctr2 = cellCentered_2ndOrder[index];
                (*this)(avg, ctr4, ctr2);
            }
        }
    };
}
#endif