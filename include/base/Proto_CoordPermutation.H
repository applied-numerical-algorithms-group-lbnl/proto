#pragma once

#ifndef _PROTO_COORD_PERMUTATION_
#define _PROTO_COORD_PERMUTATION_

#include "Proto_Point.H"
#include "Proto_Box.H"
#include "Proto_Array.H"

namespace Proto
{
    class CoordPermutation
    {
        public:

        /// Trivial Permutation
        inline static CoordPermutation identity();
        /// Clockwise (-90* rotation) in the XY plane
        inline static CoordPermutation cw(int a_coord = 2);
        /// Counter-Clockwise (+90* rotation) in the XY plane
        inline static CoordPermutation ccw(int a_coord = 2);
        /// Invert all coordinates
        inline static CoordPermutation reverse();
        /// Rotate the plane containing two coordinates such that c0 -> c1
        /**
         * Note that swapping the order of the inputs reverses the direction or rotation.
        */
        inline static CoordPermutation rotatePlane(int a_c0, int a_c1);

        /// Default / Identity constructor
        inline CoordPermutation();
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
        inline CoordPermutation(std::initializer_list<Array<int,  3>> a_args);
        inline CoordPermutation(std::initializer_list<Array<Point,2>> a_args);
        inline CoordPermutation(std::vector<Array<int, 3>>& a_args);
        inline void define(std::vector<Array<int, 3>>& a_args);
        inline void defineMatrix(Array<Array<int, DIM>, DIM> a_matrix);
        /// Direct Matrix Construction
        //CoordPermutation(Array<Array<int, DIM>, DIM> a_matrix);

        /// Permute the coordinates of a Point
        inline Point operator()(Point a_pt) const;
        
        /// Compute and return the inverse Permutation
        inline CoordPermutation inverse() const; 

        /// Rotate Cell
        /**
            Given a Point and it's containing Box in unrotated coordinates,
            compute the Point associated with the same finite volume cell in rotated coordinates.
            Note that because a finite-volume cell is designated by it's low corner, this function
            is not the same as rotatePoint.
        */
        inline Point rotateCell(Point a_point, const Box& a_srcBox, const Box& a_dstBox) const; 
        
        /// Rotate Point
        /**
            Given a Point and it's containing Box in unrotated coordinates,
            compute the Point associated with the same Point (vertex) in rotated coordinates.

        */
        inline Point rotatePoint(Point a_point, const Box& a_srcBox, const Box& a_dstBox) const; 
        
        template< typename T, unsigned int C=1, unsigned char D=1, unsigned char E=1>
        inline void rotateBuffer(T* srcData, T* dstData,
                const Box& srcBox, const Box& dstBox) const;

        inline Array<Array<int, DIM>, DIM> matrix() const {return m_matrix; }
        
        inline bool operator==(const CoordPermutation& a_rhs) const;

        inline bool isIdentity() const;

        inline void print() const;

        inline void printMatrix() const;


        private:
        inline bool isValid() const; // checks if the permutation is valid

        Array<Array<int, DIM>, DIM> m_matrix; // discrete rotation matrix
    };

    inline CoordPermutation operator*(const CoordPermutation& a_A, const CoordPermutation& a_B);
#include "implem/Proto_CoordPermutationImplem.H"
} // end namespace Proto
#endif //end include guard
