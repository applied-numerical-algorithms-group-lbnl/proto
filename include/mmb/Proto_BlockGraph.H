#pragma once

#include "Proto.H"
#include <map>
#include <array>

namespace Proto
{
    class CoordPermutation
    {
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

        Point forward(const Point& a_point);
        Point inverse(const Point& a_point);

        private:
        std::array<std::array<int, DIM>, DIM>; // discrete rotation matrix
    };

    class Block
    {
        public:
        
        Block(unsigned int a_id, unsigned int a_blockSize);
        
        const ProblemDomain& domain() const;
        Box domainBox() const;
        unsigned int id() const;

        private:

        unsigned int m_id;
        ProblemDomain m_domain;
    };

    template<unsigned int CODIM>
    class BlockBoundary
    {
        public:
        BlockBoundary();
        void addBlock(Block& a_block);
        
        private:
        std::vector<Block*> m_blocks;
    };

    template<>
    class BlockBoundary<1>
    {
        public:
        BlockBoundary(Block& a_start, Block& a_end);

        private:
        Block& m_start;
        Block& m_end;
        CoordPermutation m_permute;
    };
    
    class BlockGraph
    {
        public:

        BlockGraph(unsigned int a_numBlocks, unsigned int a_blockSize);
        
        bool addBoundary(
                unsigned int a_b0,              // starting block
                unsigned int a_b1,              // ending block
                unsigned int a_coord,           // coordinate direction
                Side::LoHiSide a_side,          // coordinate "sign"
                CoordPermutation a_rotation);   // orientation of b1 w.r.t. b0

        private:
        std::unordered_map<unsigned int, Block> m_blocks;
        std::unordered_map<Point, BlockBoundary<1>> m_boundary_1; //codim 1
        std::unordered_map<Point, BlockBoundary<2>> m_boundary_2; //codim 2
#if DIM > 2
        std::unordered_map<Point, BlockBoundary<3>> m_boundary_3; //codim 3
#endif
    };
}
