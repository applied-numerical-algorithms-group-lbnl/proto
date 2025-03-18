
#pragma once
#ifndef _PROTO_NAME_OF_YOUR_OP_
#define _PROTO_NAME_OF_YOUR_OP_

#include "Proto.H"

namespace Proto
{
    template <class T>
    class NameOfYourOp : public CompositeOp
    {

    private:
        // Stencils, intermediates, etc.

    public:
        inline NameOfYourOp()
        {
        }

        // Mandatory
        inline std::vector<Box> domains(const Box &rangeBox) const
        {
            std::vector<Box> domainBoxes;
            // compute domains of each BoxData argument of operator()
            return domainBoxes;
        }

        // Mandatory
        inline Box range(const std::vector<Box> &domainBoxes) const
        {
            // compute the range from the domains of all inputs of operator()
           return domainBoxes[0];
        }

        // Mandatory
        template <unsigned int C, MemType MEM, unsigned int D = 1, unsigned int E = 1>
        inline void
        operator()(
            BoxData<T, C, MEM, D, E> &output,
            const BoxData<T, C, MEM, D, E> &input)
        {
            PR_TIME("NameOfYourOp::operator()");
            // implement operator logic
        }

        // Optional, but easy to implement
        template <unsigned int C, MemType MEM, unsigned int D = 1, unsigned int E = 1>
        inline BoxData<T, C, MEM, D, E>
        operator()(
            const BoxData<T, C, MEM, D, E> &input)
        {
            Box rangeBox = range(getDomains(input));
            BoxData<T, C, MEM, D, E> output(rangeBox);
            (*this)(output, input);
            return output;
        }

        // Optional
        template <unsigned int C, MemType MEM, Centering CTR>
        inline void
        operator()(
            LevelBoxData<T, C, MEM, CTR> &output,
            const LevelBoxData<T, C, MEM, CTR> &input)
        {
            PROTO_ASSERT(input.layout().compatible(output.layout()),
                         "incompatible layouts");

            auto &layout = input.layout();
            for (auto index : layout)
            {
                auto &out = output[index];
                const auto &in= input[index];
                (*this)(output, input);
            }
        }
    };
}
#endif