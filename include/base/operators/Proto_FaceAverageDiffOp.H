
#pragma once
#ifndef _PROTO_FACE_DIFF_OP_
#define _PROTO_FACE_DIFF_OP_

#include "Proto.H"

namespace Proto
{
    template <class T>
    class FaceAverageDiffOp : public CompositeOp<1>
    {

    private:
        int _normalDir;
        int _order;
        Array<Stencil<T>, DIM> _faceDerivs;

    public:
        inline FaceAverageDiffOp(int normalDir, int order)
        {
            _normalDir = normalDir;
            _order = order;
            for (int dd = 0; dd < DIM; dd++)
            {
                if (dd == normalDir)
                {
                    _faceDerivs[dd] = Stencil<T>::DiffCellToFace(normalDir, Side::Lo, order);
                }
                else
                {
                    _faceDerivs[dd] =           
                        Stencil<T>::CellToFace(dd, Side::Lo, order)*
                        Stencil<T>::FluxDivergence(dd)*
                        Stencil<T>::CellToFace(normalDir, Side::Lo,order);
                }
            }
        }

        inline std::array<Box, 1> domains(const Box &rangeBox) const
        {
            std::array<Box, 1> domainBoxes;
            Box b;
            for (int dd = 0; dd < DIM; dd++)
            {
                b += _faceDerivs[dd].domain(rangeBox);
            }
            domainBoxes[0] = b;
            return domainBoxes;
        }

        inline Box range(const std::array<Box, 1> &domainBoxes) const
        {
            Box b = _faceDerivs[0].range(domainBoxes[0]);
            for (int dd = 1; dd < DIM; dd++)
            {
                b &= _faceDerivs[dd].range(domainBoxes[0]);
            }
            return b;
        }

        template <unsigned int C, MemType MEM>
        inline void
        operator()(
            BoxData<T, C, MEM, DIM> &output,
            const BoxData<T, C, MEM> &input)
        {
            PR_TIME("FaceDiff::operator()");
            for (int dir = 0; dir < DIM; dir++)
            {
                BoxData<T,C,MEM,1> grad_i = plane(output, dir);
                grad_i |= _faceDerivs[dir](input);
            }
        }

        template <unsigned int C, MemType MEM>
        inline BoxData<T, C, MEM, DIM>
        operator()(
            const BoxData<T, C, MEM> &input)
        {
            auto domains = getDomains(input);
            Box rangeBox = range(domains);
            BoxData<T, C, MEM, DIM> output(rangeBox);
            (*this)(output, input);
            return output;
        }

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