
#pragma once
#ifndef _PROTO_CONVOLVE_EDGE_OP_
#define _PROTO_CONVOLVE_EDGE_OP_

namespace Proto
{
    /// @brief Convolves data which is centered along a specified direction. The result will be averaged along that direction to 4th order.
    /// @tparam T 
    template <class T>
    class ConvolveEdgeOp : public CompositeOp<2>
    {

    private:
        Stencil<T> m_L;

    public:
        inline ConvolveEdgeOp(int dir)
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
            BoxData<T, C, MEM, D, E> &edgeAveraged_4thOrder,
            const BoxData<T, C, MEM, D, E> &edgeCentered_4thOrder,
            const BoxData<T, C, MEM, D, E> &edgeCentered_2ndOrder)
        {
            PR_TIME("ConvolveEdgeOp::operator()");
            edgeAveraged_4thOrder |= m_L(edgeCentered_2ndOrder, 1.0 / 24.0);
            edgeAveraged_4thOrder += edgeCentered_4thOrder;
        }

        template <unsigned int C, MemType MEM, unsigned int D = 1, unsigned int E = 1>
        inline BoxData<T, C, MEM, D, E>
        operator()(
            const BoxData<T, C, MEM, D, E> &edgeCentered_4thOrder,
            const BoxData<T, C, MEM, D, E> &edgeCentered_2ndOrder)
        {
            Box rangeBox = range(getDomains(edgeCentered_4thOrder, edgeCentered_2ndOrder));
            BoxData<T, C, MEM, D, E> edgeAveraged(rangeBox);
            (*this)(edgeAveraged, edgeCentered_4thOrder, edgeCentered_2ndOrder);
            return edgeAveraged;
        }

        template <unsigned int C, MemType MEM, Centering CTR>
        inline void
        operator()(
            LevelBoxData<T, C, MEM, CTR> &edgeAveraged_4thOrder,
            const LevelBoxData<T, C, MEM, CTR> &edgeCentered_4thOrder,
            const LevelBoxData<T, C, MEM, CTR> &edgeCentered_2ndOrder)
        {
            PROTO_ASSERT(edgeAveraged_4thOrder.layout().compatible(edgeCentered_4thOrder.layout()),
                         "incompatible layouts");
            PROTO_ASSERT(edgeAveraged_4thOrder.layout().compatible(edgeCentered_2ndOrder.layout()),
                         "incompatible layouts");

            auto &layout = edgeAveraged_4thOrder.layout();
            for (auto index : layout)
            {
                auto &avg = edgeAveraged_4thOrder[index];
                const auto &ctr4 = edgeCentered_4thOrder[index];
                const auto &ctr2 = edgeCentered_2ndOrder[index];
                (*this)(avg, ctr4, ctr2);
            }
        }
    };
}
#endif