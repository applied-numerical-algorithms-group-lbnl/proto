
#pragma once
#ifndef _PROTO_DECONVOLVE_FACE_OP_
#define _PROTO_DECONVOLVE_FACE_OP_

namespace Proto
{
    template <class T>
    class DeconvolveFaceOp : public CompositeOp<2>
    {

    private:
        Stencil<T> m_Lf;

    public:
        inline DeconvolveFaceOp(int dir)
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
            PROTO_ASSERT(domainBoxes.size() == 2,
                         "DeconvolveFaceOp::range | Error: DeconvolveFaceOp requires two input domains for the 4th order and 2nd order cell centered data respectively");
            return (m_Lf.range(domainBoxes[1]) & domainBoxes[0]);
        }

        template <unsigned int C, MemType MEM, unsigned int D = 1, unsigned int E = 1>
        inline void
        operator()(
            BoxData<T, C, MEM, D, E> &cellCentered_4thOrder,
            const BoxData<T, C, MEM, D, E> &cellAveraged_4thOrder,
            const BoxData<T, C, MEM, D, E> &cellAveraged_2ndOrder)
        {
            PR_TIME("DeconvolveFaceOp::operator()");
            cellCentered_4thOrder |= m_Lf(cellAveraged_2ndOrder, -1.0 / 24.0);
            cellCentered_4thOrder += cellAveraged_4thOrder;
        }

        template <unsigned int C, MemType MEM, unsigned int D = 1, unsigned int E = 1>
        inline BoxData<T, C, MEM, D, E>
        operator()(
            const BoxData<T, C, MEM, D, E> &cellAveraged_4thOrder,
            const BoxData<T, C, MEM, D, E> &cellAveraged_2ndOrder)
        {
            Box rangeBox = range(getDomains(cellAveraged_4thOrder, cellAveraged_2ndOrder));
            BoxData<T, C, MEM, D, E> cellCentered(rangeBox);
            (*this)(cellCentered, cellAveraged_4thOrder, cellAveraged_2ndOrder);
            return cellCentered;
        }

        template <unsigned int C, MemType MEM, Centering CTR>
        inline void
        operator()(
            LevelBoxData<T, C, MEM, CTR> &cellCentered_4thOrder,
            const LevelBoxData<T, C, MEM, CTR> &cellAveraged_4thOrder,
            const LevelBoxData<T, C, MEM, CTR> &cellAveraged_2ndOrder)
        {
            PROTO_ASSERT(cellCentered_4thOrder.layout().compatible(cellAveraged_4thOrder.layout()),
                         "incompatible layouts");
            PROTO_ASSERT(cellCentered_4thOrder.layout().compatible(cellAveraged_2ndOrder.layout()),
                         "incompatible layouts");

            auto &layout = cellCentered_4thOrder.layout();
            for (auto index : layout)
            {
                auto &ctr = cellCentered_4thOrder[index];
                const auto &avg4 = cellAveraged_4thOrder[index];
                const auto &avg2 = cellAveraged_2ndOrder[index];
                (*this)(ctr, avg4, avg2);
            }
        }
    };
}
#endif