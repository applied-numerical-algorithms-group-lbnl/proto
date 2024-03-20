#pragma once

#ifndef _PROTO_MB_GRAPH_
#define _PROTO_MB_GRAPH_

//#include "Proto.H"

namespace Proto
{
    typedef int BlockIndex;

    /// Codimension Directions
    /** Returns a list of points with the desired codimension in space DIM */
    inline std::vector<Point> codimDirs(unsigned int a_codim);

    /// Multiblock Graph Arc
    /** Represents a directional arc between blocks in an MBGraph */
    struct MBGraphArc
    {
        inline MBGraphArc(
                BlockIndex a_srcBlock,
                BlockIndex a_dstBlock,
                Point a_srcToDst,
                Point a_dstToSrc)
        {
            srcBlock = a_srcBlock;
            dstBlock = a_dstBlock;
            srcToDst = a_srcToDst;
            dstToSrc = a_dstToSrc;
        }

        inline ~MBGraphArc(){}

        inline void print() {
            std::cout << "srcBlock: " << srcBlock;
            std::cout << " | dstBlock: " << dstBlock;
            std::cout << " | srcToDst: " << srcToDst;
            std::cout << " | dstToSrc: " << dstToSrc << std::endl;
        }

        BlockIndex srcBlock; ///< The block index at the arc's origin
        BlockIndex dstBlock; ///< The block index at the arc's end
        Point srcToDst;        ///< Direction of dstBlock relative to srcBlock in src coords
        Point dstToSrc;        ///< Direction of srcBlock relative to dstBlock in dst coords
        CoordPermutation R;    ///< Coordinate transformation from src to dst coords
    };

    /// Multiblock Graph Node
    /** Represents a block in an MBGraph. Contains MBGraphArc objects*/
    struct MBGraphNode
    {
        inline MBGraphNode(BlockIndex a_index);
    
        /// Add Arc
        inline void addArc(BlockIndex a_block, Point a_toDst, Point a_fromDst);
        
        /// Get Coordinate Rotation
        /** Retrieve the coordinate rotation between this block and a specified block
         *  In a specified direction (in this block's coordinate system)
         *
         *  \param a_dir    Direction of desired output in this block's coordinates
         *  \param a_block  The specified neighboring block 
         */
        inline CoordPermutation& rotation(Point a_dir, BlockIndex a_block);

        BlockIndex index; ///< The block associated with this node
        mutable std::map<Point, std::vector<MBGraphArc>> arcs; ///< List of arcs
    };
  
    class MBProblemDomain;

    /// Multiblock Graph
    /** Codifies the underlying graph structure of a Multiblock domain. Contains information
     *  about how blocks are connected to each other. Does not contain information related to
     *  block sizes (see MBProblemDomain).
     *
     *  A MBGraph is built in stages. First a graph is initialized with the desired number of
     *  blocks. Then, codimension-1 boundaries are defined using the defineBoundary function.
     *  Boundaries of codimension-2+ are created automatically.
     *
     *  7-28-2023: Currently MBGraph does not correctly generate coordinate rotations associated
     *  with codimension 3+ boundaries. -CLG
     */
    class MBGraph
    {
        friend class MBProblemDomain;
        public:

        /// Constructor
        /** Creates an empty graph with a specified number of nodes and no connections */
        inline MBGraph(BlockIndex a_numBlocks);

        /// Define Boundary
        /** Adds a codimension-1 boundary to the graph. This function creates both directional
         *  arcs that are implied by the boundary (i.e. there is no need to call this function
         *  twice for each pair of connected blocks).
         *
         *  Implicitly creates any codimension-2+ boundaries. This function will result
         *  in an error if a codimension-2+ boundary or invalid geometry is specified.
         * 
         *  \param a_srcBlock   Block A
         *  \param a_dstBlock   Block B
         *  \param a_dir        Codimension-1 direction from A to B in the coords of block A
         *  \param a_rotation   Coordinate permutation from A to B
         */
        inline void defineBoundary(
                BlockIndex        a_srcBlock,
                BlockIndex        a_dstBlock,
                Point               a_dir,
                CoordPermutation&   a_rotation);

        /// Define Boundary (overload)
        inline void defineBoundary(
                BlockIndex        a_srcBlock,
                BlockIndex        a_dstBlock,
                unsigned int      a_dir,
                Side::LoHiSide      a_side,
                CoordPermutation&   a_rotation);

        /// Post Process Rotations
        /** This function is a postprocessing operation which corrects codimension-2 coordinate
         *  rotations. The function must be called after all desired boundaries are defined. 
         */
        inline void fixRotations();

        /// Size (Number of Blocks)
        inline unsigned int size() const {return m_blocks.size(); }

        /// Get Number of Boundaries
        /** If the specified boundary is codimension-1, the output will be at most 1. */
        inline unsigned int numBoundaries(BlockIndex a_blockIndex, Point& a_dir) const;
        
        /// Get Boundary Arcs
        /** Get a list of MBGraphArc objects associated with a boundary. If the boundary
         *  does not exist, an empty vector is returned. Not intended for public use.
         */
        inline std::vector<MBGraphArc> boundaries(
                BlockIndex a_srcBlock,
                Point        a_dir) const;

        /// Full Connectivity
        /** Returns the full list of directions of all codimensionalities that define a
         *  boundary between the two specified blocks. The returned Point set is in the
         *  coordinates of the source block.
         *  */
        inline std::vector<Point> fullConnectivity(
                BlockIndex a_srcBlock,
                BlockIndex a_dstBlock) const;
        
        /// Connectivity
        /** Returns the direction associated with the boundary of minimal codimension between
         *  the two specified blocks. The Point returned is in the source block coordinates.
         */
        inline Point connectivity(BlockIndex a_srcBlock, BlockIndex a_dstBlock) const;

        /// Reverse Arc
        /** Get direction from dst to src along the graph arc from a_dir. The input Point
         *  is in srcBlock's coordinates and the returned Point is in dstBlock's coordinates.
         *  The input direction must correspond to a boundary between the two blocks.
         */
        inline Point reverseArc(
                BlockIndex a_srcBlock,
                BlockIndex a_dstBlock,
                Point a_dir) const;
        
        /// Reverse Direction
        /** Mirrors a direction across a boundary accounting for coordinate transformations.
         *  The input direction is in srcBlock's coordinates and the output direction is in
         *  dstBlock's coordinates.
         */
        inline Point reverseDir(
                BlockIndex a_srcBlock,
                BlockIndex a_dstBlock,
                Point a_dir) const;

        /// Get Rotation
        /** Obtain the coordinate rotation from a specified block to an adjacent block which is
         *  connected by a specified direction. If the combination of inputs does not correspond
         *  to a valid block boundary, this function will yield an error.
         */
        inline CoordPermutation& rotation(
                BlockIndex a_srcBlock,
                Point a_dir,
                BlockIndex a_dstBlock) const;
        
        /// Get Rotation
        /** Obtain the coordinate rotation from a specified block to an adjacent block which is
         *  connected by a codimension-1 boundary. If the specified direction is not
         *  not codimension-1, this function will yield an error.
         */
        inline CoordPermutation& rotation(
                BlockIndex a_block,
                Point a_dir) const;

        /// Get Adjacent Block
        /** Obtain the block adjacent to a specified block in a specified codimension-1
         *  direction */
        inline BlockIndex adjacent(
                BlockIndex    a_srcBlock,
                unsigned int    a_dir,
                Side::LoHiSide  a_side) const;
        
        /// Get Adjacent Block (Overload)
        /** This function is the same as the other version of adjacent. This function will fail
         *  if a codimension-2+ point is supplied. */
        inline BlockIndex adjacent(
                BlockIndex    a_srcBlock,
                Point           a_dir) const;
        

        /// Equality Operator
        inline bool operator==(const MBGraph& a_rhs) const;

        /// Detect Block Boundary
        /** Returns true if the specified inputs correspond to a boundary between two blocks. */
        inline bool isBlockBoundary(BlockIndex a_srcBlock, Point a_dir) const;
      
        /// Return list of codimension-2 triple points
        inline std::vector<Point> triplePoints(
                BlockIndex a_srcBlock,
                Point a_dir) const;

        /// Detect Triple Point
        /** Returns true if the specified block boundary corresponds to a "triple point".
         *  In this context, a "triple point" is a codimension-2 block boundary where exactly
         *  three blocks meet. This function is used to infer corrections that must be made 
         *  when building interpolation operators in such regions. */
        inline bool isTriplePoint(
                BlockIndex  a_srcBlock,
                Point        a_dir) const;

        inline bool isDomainBoundary(
                BlockIndex  a_srcBlock,
                Point       a_dir) const;

        /// Print
        inline void print() const;
        
        /// Close Graph
        /** Executes postprocessing steps once all boundaries have been defined. This function
         *  Also checks for various topological errors in the graph construction. */
        inline void close();

        private:
        bool m_closed;
        std::vector<std::shared_ptr<MBGraphNode>> m_blocks;
        
        
    };
#include "implem/Proto_MBGraphImplem.H"
} // end namespace Proto
#endif // end include guard
