
#pragma once
#ifndef _PROTO_MAPPED_FACE_AVERAGE_DIFF_OP_
#define _PROTO_MAPPED_FACE_AVERAGE_DIFF_OP_

#include "Proto.H"

namespace Proto
{
    template <class T>
    class MappedFaceAverageDiffOp : public CompositeOp
    {

    private:
        int _normDir;
        std::unique_ptr<FaceAverageDiffOp<T>> _diff2Op;
        std::unique_ptr<FaceAverageDiffOp<T>> _diff4Op;
        std::unique_ptr<FaceAverageMatrixProductOp<T>> _matProdOp;
        std::unique_ptr<FaceAverageTensorQuotientOp<T>> _quotientOp;

    public:
        inline MappedFaceAverageDiffOp(int normDir) : _normDir(normDir)
        {
            _diff2Op = std::make_unique<FaceAverageDiffOp<T>>(normDir, 2);
            _diff4Op = std::make_unique<FaceAverageDiffOp<T>>(normDir, 4);
            _matProdOp = std::make_unique<FaceAverageMatrixProductOp<T>>(normDir);
            _quotientOp = std::make_unique<FaceAverageTensorQuotientOp<T>>(normDir);
        }

        inline std::vector<Box> domains(const Box &rangeBox) const
        {
            auto [jgradphi4, J4, jgradphi2, J2] = 
                tupleFromVector(_quotientOp->domains(rangeBox));
            auto [gradxphi4, N4, gradxphi2, N2] = 
                tupleFromVector(_matProdOp->domains(jgradphi4));
            auto [phi4] = tupleFromVector(_diff4Op->domains(gradxphi4));
            auto gradxphi2 = jgradphi2;
            N2 += jgradphi2;
            auto [phi2] = tupleFromVector(_diff2Op->domains(gradxphi2));

            std::vector<Box> output;
            output.push_back(phi4);
            output.push_back(phi2);
            output.push_back(N4);
            output.push_back(N2);
            output.push_back(J4);
            output.push_back(J2);
            return output;
        }

        inline Box range(const std::vector<Box> &domainBoxes) const
        {
            auto [phi4, phi2, N4, N2, J4, J2] = tupleFromVector(domainBoxes);
            auto G4 = _diff4Op->range(phi4);
            auto G2 = _diff2Op->range(phi2);
            auto JG4 = _matProdOp->range(G4, N4, G2, N2);
            auto JG2 = (G2 & N2);
            return _quotientOp->range(JG4, J4, JG2, J2);
        }
        
        template <unsigned int C, MemType MEM>
        inline void
        operator()(
            BoxData<T, C, MEM, DIM>& diff_4thOrder,
            const BoxData<T, C, MEM, 1> &phiFace_4thOrder,
            const BoxData<T, C, MEM, 1> &phiFace_2ndOrder,
            const BoxData<T, DIM, MEM, DIM> &NFace_4thOrder,
            const BoxData<T, DIM, MEM, DIM> &NFace_2ndOrder,
            const BoxData<T, 1, MEM, 1> &JFace_4thOrder,
            const BoxData<T, 1, MEM, 1> &JFace_2ndOrder)
        {
            PR_TIME("MappedFaceAverageDiffOp::operator()");
            auto G4 = (*_diff4Op)(phiFace_4thOrder);
            auto G2 = (*_diff2Op)(phiFace_2ndOrder);
            auto JG4 = (*_matProdOp)(G4, NFace_4thOrder, G2, NFace_2ndOrder);
            auto JG2 = matrixProduct(G2, NFace_2ndOrder);
            (*_quotientOp)(diff_4thOrder, JG4, JFace_4thOrder, JG2, JFace_2ndOrder);
        }

        template <unsigned int C, MemType MEM, unsigned int D = 1, unsigned int E = 1>
        inline BoxData<T, C, MEM, DIM>
        operator()(
            const BoxData<T, C, MEM, 1> &phiFace_4thOrder,
            const BoxData<T, C, MEM, 1> &phiFace_2ndOrder,
            const BoxData<T, DIM, MEM, DIM> &NFace_4thOrder,
            const BoxData<T, DIM, MEM, DIM> &NFace_2ndOrder,
            const BoxData<T, 1, MEM, 1> &JFace_4thOrder,
            const BoxData<T, 1, MEM, 1> &JFace_2ndOrder)
        {
            Box rangeBox = range(getDomains(
                phiFace_4thOrder,
                phiFace_2ndOrder,
                NFace_4thOrder,
                NFace_2ndOrder,
                JFace_4thOrder,
                JFace_2ndOrder
            ));
            BoxData<T, C, MEM, D, E> output(rangeBox);
            (*this)(output,
                phiFace_4thOrder, phiFace_2ndOrder,
                NFace_4thOrder, NFace_2ndOrder,
                JFace_4thOrder, JFace_2ndOrder);
            return output;
        }
    };
}
#endif