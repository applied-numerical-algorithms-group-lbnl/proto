
#pragma once
#ifndef _PROTO_DECONVOLVE_EDGE_OP_
#define _PROTO_DECONVOLVE_EDGE_OP_

namespace Proto
{
    template <class T>
    class DeconvolveEdgeOp : public CompositeOp<2>
    {

    private:
        Stencil<T> m_L;

    public:
        inline DeconvolveEdgeOp(int dir)
        {
            m_L = (+1.0)*Shift::Basis(dir,-1)
                + (-2.0)*Shift::Zeros()
                + (+1.0)*Shift::Basis(dir,1);
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
            return (m_L.range(domainBoxes[1]) & domainBoxes[0]);
        }

        template <unsigned int C, MemType MEM, unsigned int D = 1, unsigned int E = 1>
        inline void
        operator()(
            BoxData<T, C, MEM, D, E> &edgeCentered_4thOrder,
            const BoxData<T, C, MEM, D, E> &edgeAveraged_4thOrder,
            const BoxData<T, C, MEM, D, E> &edgeAveraged_2ndOrder)
        {
            PR_TIME("DeconvolveEdgeOp::operator()");
            edgeCentered_4thOrder |= m_L(edgeAveraged_2ndOrder, -1.0 / 24.0);
            edgeCentered_4thOrder += edgeAveraged_4thOrder;
        }

        template <unsigned int C, MemType MEM, unsigned int D = 1, unsigned int E = 1>
        inline BoxData<T, C, MEM, D, E>
        operator()(
            const BoxData<T, C, MEM, D, E> &edgeAveraged_4thOrder,
            const BoxData<T, C, MEM, D, E> &edgeAveraged_2ndOrder)
        {
            Box rangeBox = range(getDomains(edgeAveraged_4thOrder, edgeAveraged_2ndOrder));
            BoxData<T, C, MEM, D, E> edgeCentered(rangeBox);
            (*this)(edgeCentered, edgeAveraged_4thOrder, edgeAveraged_2ndOrder);
            return edgeCentered;
        }

        template <unsigned int C, MemType MEM, Centering CTR>
        inline void
        operator()(
            LevelBoxData<T, C, MEM, CTR> &edgeCentered_4thOrder,
            const LevelBoxData<T, C, MEM, CTR> &edgeAveraged_4thOrder,
            const LevelBoxData<T, C, MEM, CTR> &edgeAveraged_2ndOrder)
        {
            PROTO_ASSERT(edgeCentered_4thOrder.layout().compatible(edgeAveraged_4thOrder.layout()),
                         "incompatible layouts");
            PROTO_ASSERT(edgeCentered_4thOrder.layout().compatible(edgeAveraged_2ndOrder.layout()),
                         "incompatible layouts");

            auto &layout = edgeCentered_4thOrder.layout();
            for (auto index : layout)
            {
                auto &ctr = edgeCentered_4thOrder[index];
                const auto &avg4 = edgeAveraged_4thOrder[index];
                const auto &avg2 = edgeAveraged_2ndOrder[index];
                (*this)(ctr, avg4, avg2);
            }
        }
    };
}
#endif