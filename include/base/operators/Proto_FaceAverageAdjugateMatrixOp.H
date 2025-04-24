
#pragma once
#ifndef _PROTO_FACE_AVERAGE_ADJUGATE_MATRIX_OP_
#define _PROTO_FACE_AVERAGE_ADJUGATE_MATRIX_OP_

#include "Proto.H"

namespace Proto
{
    template <class T>
    class FaceAverageAdjugateMatrixOp : public CompositeOp<1>
    {

    private:
        Array<std::unique_ptr<FaceAverageSurfaceTransformOp<T,DIM>>, DIM> _columnOps;
        Array<Stencil<T>,DIM> _avgOps;
        int _normDir;

    public:
        inline FaceAverageAdjugateMatrixOp(int normDir)
        {
            _normDir = normDir;
            for (int dd = 0; dd < DIM; dd++)
            {
                _columnOps[dd] = std::make_unique<FaceAverageSurfaceTransformOp<T,DIM>>(dd);
                if (dd != normDir)
                {
                    _avgOps[dd] = Stencil<T>::faceToCell(dd)*Stencil<T>::CellToFace(normDir);
                }
            }
        }

        inline std::array<Box, 1> domains(const Box &rangeBox) const
        {
            std::array<Box, 1> domainBoxes;
            domainBoxes[0] = rangeBox;
            for (int jj = 0; jj < DIM; jj++)
            {
                Box bj = rangeBox;
                if (jj != _normDir)
                {
                    bj = _avgOps[jj].domain(rangeBox);
                }
                domainBoxes[0] += (*_columnOps[jj]).domains(bj)[0];
            }
            return domainBoxes;
        }

        inline Box range(const std::array<Box, 1> &domainBoxes) const
        {
            Box domain = domainBoxes[0];
            Box range = domain;
            for (int jj = 0; jj < DIM; jj++)
            {
                Box adj_j = _columnOps[jj]->range({domain});
                if (jj != _normDir)
                {
                    adj_j &= _avgOps[jj].range(adj_j);
                }
                range &= adj_j;
            }

            return range;
        }

        template <MemType MEM>
        inline void
        operator()(
            FluxBoxData<T, DIM, MEM> &adjugateComponents,
            const BoxData<T, DIM, MEM> &nodeCoords)
        {
            PR_TIME("FaceAverageAdjugateMatrixOp::operator(FluxBoxData, ...)");
            for (int dd = 0; dd < DIM; dd++)
            {
                (*_columnOps[dd])(adjugateComponents[dd], nodeCoords);
            }
        }

        template <MemType MEM>
        inline void
        operator()(
            BoxData<T, DIM, MEM, DIM> &adjugate,
            const BoxData<T, DIM, MEM> &nodeCoords)
        {
            PR_TIME("FaceAverageAdjugateMatrixOp::operator(BoxData, ...)");
            for (int jj = 0; jj < DIM; jj++) //column
            {
                BoxData<T,DIM,MEM,1> adj_j = plane(adjugate, jj);
                if (jj != _normDir)
                {
                    auto tmp = (*_columnOps[jj])(nodeCoords);
                    adj_j |= _avgOps[jj](tmp);
                } else {
                    (*_columnOps[jj])(adj_j, nodeCoords);
                }
            }
        }

        template <MemType MEM>
        inline BoxData<T, DIM, MEM, DIM>
        operator()(
            const BoxData<T, DIM, MEM> &nodeCoords)
        {
            Box rangeBox = range(getDomains(nodeCoords));
            BoxData<T, DIM, MEM, DIM> output(rangeBox);
            (*this)(output, nodeCoords);
            return output;
        }
    };
}
#endif