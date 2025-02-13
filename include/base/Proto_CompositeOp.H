#pragma once
#ifndef _PROTO_COMPOSITE_OP_
#define _PROTO_COMPOSITE_OP_

#include "Proto.H"

namespace Proto 
{
    class CompositeOp
    {
    public:
        inline virtual std::vector<Box> domains(const Box &rangeBox) const = 0;
        inline virtual Box range(const std::vector<Box> &domainBoxes) const = 0;
        template<typename... Args>
        inline Box getRange(Box first, Args... boxes)
        {
            static_assert((std::is_same_v<Box, Args>&& ...), "All arguments must be of type Box");
            std::vector<Box> domains;
            domains.push_back(first);
            (domains.push_back(boxes), ...);
            return range(domains);
        }

    protected:
        
        template<typename... Srcs>
        inline std::vector<Box> getDomains(Srcs... sources)
        {
            std::vector<Box> domains;
            (domains.push_back(sources.box()), ...); //C++17 fold expression
            return domains;
        }
    };
}
#include "operators/Proto_ConvolveOp.H"
#include "operators/Proto_DeconvolveOp.H"
#include "operators/Proto_ConvolveFaceOp.H"
#include "operators/Proto_DeconvolveFaceOp.H"
#endif