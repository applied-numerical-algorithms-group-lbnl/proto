
#pragma once
#ifndef _PROTO_FACE_AVERAGED_SURFACE_TRANSFORM_OP_
#define _PROTO_FACE_AVERAGED_SURFACE_TRANSFORM_OP_

#include "Proto.H"

namespace Proto
{
    /// @brief Computes the value of a unit normal vector in a given direction under a specified coordinate transformation. \
    The resulting vector is unnormalized, and averaged over surfaces corresponding to the specified normal vector. \
    This quantity is equivalent to the Nth column of the adjugate of the Jacobian matrix.
    /// @tparam T 
    /// @tparam NCOMPS 
    template <class T, unsigned int NCOMPS>
    class FaceAveragedSurfaceTransformOp : public CompositeOp<1>
    {

    private:
        // Stencils, intermediates, etc.

    public:
        inline FaceAveragedSurfaceTransformOp()
        {
            // unspecialized template is not valid
            PROTO_ASSERT(false, "Error: FaceAverageSurfaceTransformOp is not implemented for NCOMPS = %u", NCOMPS);
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
    class FaceAveragedSurfaceTransformOp<T, 2> : public CompositeOp<1>
    {
        private:
        int _normDir;
        int _perpDir;
        Stencil<T> _diffOp;

        public:

        inline FaceAveragedSurfaceTransformOp(int normDir)
        {
            _normDir = normDir;
            _perpDir = (normDir + 1) % 2;
            _diffOp = Stencil<T>::FluxDivergence(_perpDir);
        }

        inline std::array<Box, 1> domains(const Box &rangeBox) const
        {
            return {_diffOp.domain(rangeBox)};
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
            PR_TIME("FaceAverageSurfaceTransformOp::operator()");

            for (int ii = 0; ii < 2; ii++)
            {
                int jj = (ii+1) % 2;
                auto Ci = slice(cofactors, ii);
                auto Xj = slice(nodeCoords, jj);
                T sign = (jj == _perpDir) ? 1.0 : -1.0;
                Ci |= _diffOp(Xj, sign);
            }
        }

        inline BoxData<T, 2>
        operator()(
            const BoxData<T, 2> &nodeCoords)
        {
            auto domains = getDomains(nodeCoords);
            Box rangeBox = range(domains);
            BoxData<T, 2> cofactors(rangeBox);
            (*this)(cofactors, nodeCoords);
            return cofactors;
        }

        template<MemType MEM>
        inline void
        operator()(
            LevelBoxData<T, 2, MEM, PR_CELL> &cofactors,
            const LevelBoxData<T, 2, MEM, PR_NODE> &nodeCoords)
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
    class FaceAveragedSurfaceTransformOp<T, 3> : public CompositeOp<1>
    {
        private:
        int _normDir;
        int _perpDir0;
        int _perpDir1;
        std::unique_ptr<EdgeAverageCrossProductOp<T>> _cross0;
        std::unique_ptr<EdgeAverageCrossProductOp<T>> _cross1;
        Stencil<T> _diff0;
        Stencil<T> _diff1;
        Stencil<T> _avg0;
        Stencil<T> _avg1;


        public:

        inline FaceAveragedSurfaceTransformOp(int normDir)
        {
            _normDir = normDir;
            _perpDir0 = (normDir + 1) % 3;
            _perpDir1 = (normDir + 2) % 3;
            _cross0 = std::make_unique<EdgeAverageCrossProductOp<T>>(_perpDir0);
            _cross1 = std::make_unique<EdgeAverageCrossProductOp<T>>(_perpDir1);
            _diff0 = Stencil<T>::FluxDivergence(_perpDir0);
            _diff1 = Stencil<T>::FluxDivergence(_perpDir1);
            _avg0 = Stencil<T>::faceToCell(_perpDir0, 4);
            _avg1 = Stencil<T>::faceToCell(_perpDir1, 4);
        }

        inline std::array<Box, 1> domains(const Box &rangeBox) const
        {
            Box c1 = _diff0.domain(rangeBox);
            Box c0 = _diff1.domain(rangeBox);
            auto [a0, d0, _a0, _d0] = _cross0->domains(c0);
            auto [a1, d1, _a1, _d1] = _cross1->domains(c1);
            a0 += _a0;
            d0 += _d0;
            a1 += _a1;
            d1 += _d1;
            Box domain;
            domain += _avg0.domain(a0);
            domain += _avg1.domain(a1);
            domain += _diff0.domain(d0);
            domain += _diff1.domain(d1);
            return {domain};
        }

        inline Box range(const std::array<Box, 1> &domainBoxes) const
        {
            Box domain = domainBoxes[0];
            Box a0 = _avg0.range(domain);
            Box a1 = _avg1.range(domain);
            Box d0 = _diff0.range(domain);
            Box d1 = _diff1.range(domain);
            Box c0 = _cross0->range({d0, a0, d0, a0});
            Box c1 = _cross1->range({d1, a1, d1, a1});
            Box range = _diff1.range(c0);
            range &= _diff0.range(c1);
            return range;
        }

        inline void
        operator()(
            BoxData<T, 3> &cofactors,
            const BoxData<T, 3> &nodeCoords)
        {
            PR_TIME("FaceAverageSurfaceTransformOp::operator()");
            BoxData<T,3> xAvg0 = _avg0(nodeCoords); 
            BoxData<T,3> xAvg1 = _avg1(nodeCoords); 
            BoxData<T,3> dx0 = _diff0(nodeCoords);  
            BoxData<T,3> dx1 = _diff1(nodeCoords);              
            BoxData<T,3> dx0_cross_x = (*_cross0)(dx0, xAvg0, dx0, xAvg0);
            BoxData<T,3> x_cross_dx1 = (*_cross1)(xAvg1, dx1, xAvg1, dx1);
            cofactors |= _diff1(dx0_cross_x, 0.5);
            cofactors += _diff0(x_cross_dx1, 0.5);
        }

        inline BoxData<T, 3>
        operator()(
            const BoxData<T, 3> &nodeCoords)
        {
            Box rangeBox = range(getDomains(nodeCoords));
            BoxData<T, 3> cofactors(rangeBox);
            (*this)(cofactors, nodeCoords);
            return cofactors;
        }

        template<MemType MEM>
        inline void
        operator()(
            LevelBoxData<T, 3, MEM, PR_CELL> &cofactors,
            const LevelBoxData<T, 3, MEM, PR_NODE> &nodeCoords)
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