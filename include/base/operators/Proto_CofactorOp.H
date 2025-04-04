
#pragma once
#ifndef _PROTO_NAME_OF_YOUR_OP_
#define _PROTO_NAME_OF_YOUR_OP_

#include "Proto.H"

namespace Proto
{
    template <class T, unsigned int NCOMPS>
    class CofactorOp : public CompositeOp<1>
    {

    private:
        // Stencils, intermediates, etc.

    public:
        inline CofactorOp()
        {
            // unspecialized template is not valid
            PROTO_ASSERT(false, "Error: CofactorOp is not implemented for NCOMPS = %u", NCOMPS);
        }

        inline std::array<Box, 1> domains(const Box &rangeBox) const
        {
            // should never be called
            std::array<Box, 1> domainBoxes;
            return domainBoxes;
        }

        inline Box range(const std::array<Box, 1> &domainBoxes) const
        {
            // should never be called
           return domainBoxes[0];
        }
    };

    template <class T>
    class CofactorOp<T, 2> : public CompositeOp<1>
    {
        private:
        int _normDir;
        int _perpDir;
        Stencil<T> _diffOp;

        public:

        inline CofactorOp(int normDir)
        {
            _normDir = normDir;
            _perpDir = (normDir + 1) % 2;
            _diffOp = Stencil<T>::FluxDivergence(_perpDir)
        }

        inline std::array<Box, 1> domains(const Box &rangeBox) const
        {
            return {rangeBox.grow(PR_NODE)};
        }

        inline Box range(const std::array<Box, 1> &domainBoxes) const
        {
            return _diffOp.range(domainBoxes[0]);
        }

        inline void
        operator()(
            BoxData<T, 2> &cofactors,
            const BoxData<T, 2> &nodeCoords)
        {
            PR_TIME("CofactorOp::operator()");

            int normDir = _diffDir;
            int orthDir = (_diffDir + 1) % 2;
            Box bxret = a_X.box()&a_X.box().shift(dirperp,-1);
            BoxData<T,DIM,MEM> ret(bxret);
            ret.setToZero();

            for (int comp = 0;comp < DIM;comp++)
            {
                int compperp = (comp+1)%DIM;
                auto retslice = slice(ret,comp);
                auto Xslice = slice(a_X,compperp);
                T sign = 1.0;
                if (compperp!=dirperp) sign = -1.0;
                retslice += Stencil<T>::FluxDivergence(dirperp)(Xslice,sign);
            }
            return ret;
        }

        inline BoxData<T, 2>
        operator()(
            const BoxData<T, 2> &nodeCoords)
        {
            Box rangeBox = range(getDomains(nodeCoords));
            BoxData<T, 2> cofactors(rangeBox);
            (*this)(cofactors, nodeCoords);
            return cofactors;
        }

        template<MemType MEM>
        inline void
        operator()(
            LevelBoxData<T, 2, MEM, PR_CELL> &cofactors,
            const LevelBoxData<T, C, MEM, PR_NODE> &nodeCoords)
        {
            PROTO_ASSERT(nodeCoords.layout().compatible(cofactors.layout()),
                         "incompatible layouts");

            auto &layout = nodeCoords.layout();
            for (auto index : layout)
            {
                auto &out = cofactors[index];
                const auto &in= nodeCoords[index];
                (*this)(out, in);
            }
        }
    };

    template <class T>
    class CofactorOp<T, 3> : public CompositeOp<1>
    {
        private:
        int _diffDir;

        public:

        inline CofactorOp(int diffDir) : _diffDir(diffDir)
        {
        }

        inline std::array<Box, 1> domains(const Box &rangeBox) const
        {
        }

        inline Box range(const std::array<Box, 1> &domainBoxes) const
        {
        }

        inline void
        operator()(
            BoxData<T, 2> &cofactors,
            const BoxData<T, 2> &nodeCoords)
        {
            PR_TIME("CofactorOp::operator()");
            
        }

        inline BoxData<T, 2>
        operator()(
            const BoxData<T, 2> &nodeCoords)
        {
            Box rangeBox = range(getDomains(nodeCoords));
            BoxData<T, 2> cofactors(rangeBox);
            (*this)(cofactors, nodeCoords);
            return cofactors;
        }

        template<MemType MEM>
        inline void
        operator()(
            LevelBoxData<T, 2, MEM, PR_CELL> &cofactors,
            const LevelBoxData<T, C, MEM, PR_NODE> &nodeCoords)
        {
            PROTO_ASSERT(nodeCoords.layout().compatible(cofactors.layout()),
                         "incompatible layouts");

            auto &layout = nodeCoords.layout();
            for (auto index : layout)
            {
                auto &out = cofactors[index];
                const auto &in= nodeCoords[index];
                (*this)(out, in);
            }
        }
    };

}
#endif