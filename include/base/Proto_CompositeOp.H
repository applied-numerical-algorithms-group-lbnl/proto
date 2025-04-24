#pragma once
#ifndef _PROTO_COMPOSITE_OP_
#define _PROTO_COMPOSITE_OP_

#include "Proto.H"

namespace Proto 
{
    template<typename T, std::size_t... I>
    auto tupleFromVector(const std::vector<T>& v, std::index_sequence<I...>) {
        return std::make_tuple(v[I]...);
    }

    template<typename T, std::size_t N>
    auto tupleFromVector(const std::vector<T>& v) {
        PROTO_ASSERT(v.size() < N,
            "Vector size %u is smaller than expected size %u", v.size(), N);
        return tupleFromVector<T>(v, std::make_index_sequence<N>{});
    }
    template<unsigned int NDOMAINS>
    class CompositeOp
    {
    public:
        
        inline virtual std::array<Box, NDOMAINS> domains(const Box &rangeBox) const = 0;
        
        inline virtual Box range(const std::array<Box, NDOMAINS> &domainBoxes) const = 0;
        
        template<typename... Args>
        inline Box getRange(Box first, Args... boxes)
        {
            static_assert((std::is_same_v<Box, Args>&& ...), "All arguments must be of type Box");
            std::vector<Box> domains;
            domains.push_back(first);
            (domains.push_back(boxes), ...);
            PROTO_ASSERT(domains.size() == NDOMAINS,
                "Incorrect number of domains (%u should be %u)", domains.size(), NDOMAINS);
            std::array<Box, NDOMAINS> domainArray;
            std::copy_n(domains.begin(), NDOMAINS, domainArray.begin());
            return range(domainArray);
        }
        inline Box domainUnion(const Box& rangeBox) const
        {
            Box unionBox;
            for (auto B : domains(rangeBox))
            {
                unionBox += B;
            }
            return unionBox;
        }

    protected:
        
        template<typename... Srcs>
        inline std::array<Box, NDOMAINS> getDomains(Srcs&... sources)
        {
            std::vector<Box> domains;
            (domains.push_back(sources.box()), ...); //C++17 fold expression
            PROTO_ASSERT(domains.size() == NDOMAINS,
                "incorrect number of arguments (detected %u, should be %u)",
            domains.size(), NDOMAINS);
            std::array<Box, NDOMAINS> domainArray;
            std::copy_n(domains.begin(), NDOMAINS, domainArray.begin());
            return domainArray;
        }
    };
}
#include "operators/Proto_ConvolveOp.H"
#include "operators/Proto_DeconvolveOp.H"
#include "operators/Proto_ConvolveFaceOp.H"
#include "operators/Proto_DeconvolveFaceOp.H"
#include "operators/Proto_ConvolveEdgeOp.H"
#include "operators/Proto_DeconvolveEdgeOp.H"

#include "operators/Proto_FaceAverageDiffOp.H"
#include "operators/Proto_ProductOrQuotientAverageOp.H"
#include "operators/Proto_ProductAverageOp.H"
#include "operators/Proto_CellAverageProductOp.H"
#include "operators/Proto_FaceAverageProductOp.H"
#include "operators/Proto_EdgeAverageProductOp.H"
#include "operators/Proto_QuotientAverageOp.H"
#include "operators/Proto_CellAverageQuotientOp.H"
#include "operators/Proto_FaceAverageQuotientOp.H"

#include "operators/Proto_MatrixProductOp.H"
#include "operators/Proto_FaceAverageMatrixProductOp.H"
#include "operators/Proto_CellAverageMatrixProductOp.H"
#include "operators/Proto_FaceAverageTensorQuotientOp.H"

#include "operators/Proto_CrossProductOp.H"
#include "operators/Proto_EdgeAverageCrossProductOp.H"

#include "operators/Proto_FaceAverageSurfaceTransformOp.H"
#include "operators/Proto_FaceAverageAdjugateMatrixOp.H"
#include "operators/Proto_MappedFaceAverageDiffOp.H"
#endif