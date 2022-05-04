#pragma once

#ifndef _PROTO_BLOCK_GRAPH_
#define _PROTO_BLOCK_GRAPH_

#include "Proto.H"
#include <array>

namespace Proto
{
    class CoordPermutation
    {

        /// Trivial Permutation
        static CoordPermutation identity();
        /// Clockwise (-90* rotation) in the XY plane
        static CoordPermutation cw();
        /// Counter-Clockwise (+90* rotation) in the XY plane
        static CoordPermutation ccw();
        /// 180* rotation in the XY plane
        static CoordPermutation reverse();

        /// Default / Identity constructor
        CoordPermutation(){};
        /// Constructor
        /** Syntax:
         *  <code>CoordPermutation p {{c0, c1, s}, {c0, c1, s}};</code>
         *  Which is meant to be read as "map coordinate c0 to coordinate c1 with
         *  an optional reflection s".
         *  Any coordinates omitted from the constructor are understood to be mapped
         *  through a positive identity
         *
         *  Example: 2D rotation in 3D:
         *  <code>CoordPermutation rot{{0,1,1}, {1,0,-1}};</code>
         *  coord 0 -> +coord 1
         *  coord 1 -> -coord 0
         *  coord 2 -> +coord 2
        */
        CoordPermutation(std::initializer_list<std::array<int, 3>> a_args);
        
        /// Direct construction for old school folks
        CoordPermutation(std::array<std::array<int, DIM>, DIM> a_matrix);

        /// Permute the coordinates of a Point
        Point operator()(Point a_pt) const;
        
        /// Compute and return the inverse Permutation
        CoordPermutation inverse() const; 

        private:
        bool isValid() const; // checks if the permutation is valid

        std::array<std::array<int, DIM>, DIM>; // discrete rotation matrix
    };
#include "implem/Proto_CoordPermutationImplem.H"
} // end namespace Proto
#endif //end include guard
