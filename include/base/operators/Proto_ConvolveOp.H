
#pragma once
#ifndef _PROTO_CONVOLVE_OP_
#define _PROTO_CONVOLVE_OP_

namespace Proto
{
    template <class T>
    class ConvolveOp : public CompositeOp<2>
    {

    private:
        Stencil<T> m_L;

    public:
        inline ConvolveOp()
        {
            m_L = Stencil<T>::Laplacian();
        }

        inline std::array<Box, 2> domains(const Box &rangeBox) const
        {
            std::array<Box, 2> domainBoxes;
            domainBoxes[0] = rangeBox;
            domainBoxes[1] = m_L.domain(rangeBox);
            return domainBoxes;
        }

        inline Box range(const std::array<Box, 2> &domainBoxes) const
        {
            PROTO_ASSERT(domainBoxes.size() == 2,
                         "ConvolveOp::range | Error: ConvolveOp requires two input domains for the 4th order and 2nd order cell centered data respectively");
            return m_L.range(domainBoxes[1]) & domainBoxes[0];
        }

        template <unsigned int C, MemType MEM, unsigned int D = 1, unsigned int E = 1>
        inline void
        operator()(
            BoxData<T, C, MEM, D, E> &cellAveraged_4thOrder,
            const BoxData<T, C, MEM, D, E> &cellCentered_4thOrder,
            const BoxData<T, C, MEM, D, E> &cellCentered_2ndOrder)
        {
            PR_TIME("ConvolveOp::operator()");
            Box rangeBox = range({cellCentered_4thOrder.box(), cellCentered_2ndOrder.box()});
            BoxData<T,C,MEM,D,E> tmp = m_L(cellCentered_2ndOrder, rangeBox, 1.0 / 24.0);
            tmp += cellCentered_4thOrder;
            tmp.copyTo(cellAveraged_4thOrder);
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