#pragma once

#ifndef _PROTO_MB_GRAPH_
#define _PROTO_MB_GRAPH_

namespace Proto
{
    /// @brief Defines what type is used for indexing block entities
    typedef int BlockIndex;

    /// @brief An internal structure defining an arc in the MBGraph
    struct MBGraphArc
    {
        inline MBGraphArc(
                BlockIndex srcBlock,
                BlockIndex dstBlock,
                Point srcToDst,
                Point dstToSrc)
        {
            this->srcBlock = srcBlock;
            this->dstBlock = dstBlock;
            this->srcToDst = srcToDst;
            this->dstToSrc = dstToSrc;
        }

        inline ~MBGraphArc(){}

        inline bool operator==(const MBGraphArc& rhs) const
        {
            return (srcBlock == rhs.srcBlock)
                && (dstBlock == rhs.dstBlock)
                && (srcToDst == rhs.srcToDst)
                && (dstToSrc == rhs.dstToSrc)
                && (R        == rhs.R);
        } 

        std::string str() const
        {
            std::ostringstream oss;
            oss << "MBGraphArc";
            oss << " | srcBlock: " << srcBlock;
            oss << " | dstBlock: " << dstBlock;
            oss << " | srcToDst: " << srcToDst;
            oss << " | dstToSrc: " << dstToSrc;
            return oss.str();
        }  

        inline void print(std::ostream& stream = std::cout) const { stream << str() << std::endl; }

        BlockIndex srcBlock; ///< The block index at the arc's origin
        BlockIndex dstBlock; ///< The block index at the arc's end
        Point srcToDst;        ///< Direction of dstBlock relative to srcBlock in src coords
        Point dstToSrc;        ///< Direction of srcBlock relative to dstBlock in dst coords
        CoordPermutation R;    ///< Coordinate transformation from src to dst coords
    };

    /// @brief An internal structure defining a node in the MBGraph
    struct MBGraphNode
    {
        inline MBGraphNode(BlockIndex index);
    
        /// @brief Adds an arc to this node, ingoring self and redundant arcs.
        /// The new arc represents a connection to block along direction toDst
        /// with a reverse arc direction fromDst (from dst to this). Note that
        /// toDst is in the local coordinates and fromDst is in the dst coordinates
        inline void addArc(BlockIndex block, Point toDst, Point fromDst);
        
        /// @brief Get the rotation object associated with the arc from 
        /// this along dir to block
        inline CoordPermutation& rotation(Point dir, BlockIndex block);

        /// @brief Equality established by value
        inline bool operator==(const MBGraphNode& rhs) const
        {
            return (index == rhs.index)
                && (arcs == rhs.arcs);
        }

        /// @brief Print the block associated with this node and any arcs extending from it
        /// @param stream 
        inline void print(std::ostream& stream = std::cout) const;

        BlockIndex index; ///< The block associated with this node
        mutable std::map<Point, std::vector<MBGraphArc>> arcs; /// Map from direction from this to all arcs in that direction
    };
  
    class MBProblemDomain;

    /// @brief Graph of a mapped multiblock domain. Nodes represent blocks and arcs represent
    /// boundaries between blocks. Does not contain any information about the size of blocks; 
    /// that information is provided by MBProblemDomain.
    /// 
    /// MBGraph stores connections of all codimensionalities between connected blocks.
    /// This means that if two blocks are connected
    /// by a codimension 1 boundary (sharing faces in R3) then the MBGraph also stores the 
    /// codimension 2 boundaries associated with the faces' edges, codimension 3 boundaries associated
    /// with the faces' corners, etc. up to codimension = DIM. 
    ///
    /// MBGraph currently does not completely support configurations wherein two blocks form a topological 
    /// cycle. For example, two blocks cannot have multiple codimension-1 boundaries with eachother; there
    /// must be at least one other block in the circuit for the structure to be valid.
    /// 
    /// MBGraphs usually do not need to be constructed directly, but are created implicitly when a user
    /// constructs an MBProblemDomain.
    class MBGraph
    {
        friend class MBProblemDomain;
        public:

        /// Constructor
        /** Creates an empty graph with a specified number of nodes and no connections */
        inline MBGraph(BlockIndex numBlocks);

        /// Define Boundary
        /** Adds a codimension-1 boundary to the graph. This function creates both directional
         *  arcs that are implied by the boundary (i.e. there is no need to call this function
         *  twice for each pair of connected blocks), but doing so will not corrupt the graph
         *
         *  Implicitly creates any codimension-2+ boundaries. This function will result
         *  in an error if a codimension-2+ boundary or invalid geometry is specified. 
         * 
         *  Cannot be called after this->close() has been called
         *  
         *  \param srcBlock   Block A
         *  \param dstBlock   Block B
         *  \param dir        Codimension-1 direction from A to B in the coords of block A
         *  \param rotation   Coordinate permutation from A to B
         */
        inline void defineBoundary(
                BlockIndex        srcBlock,
                BlockIndex        dstBlock,
                Point               dir,
                CoordPermutation&   rotation);

        /// Define Boundary overload using dir + side syntax instead of Point
        inline void defineBoundary(
                BlockIndex        srcBlock,
                BlockIndex        dstBlock,
                unsigned int      dir,
                Side::LoHiSide      side,
                CoordPermutation&   rotation);

        /** @brief executes post processing steps once all boundaries of this graph have been defined
         *  This function must be called before the graph can be used.
         * 
         *  Many functions of MBGraph will fail if the graph has not been closed.
         * 
         *  One of the purposes of this function is to implicitly compute higher codimensional boundaries
         *  between blocks and the coordinate permutations associated with those connections. This process
         *  is correct for codimension-2 boundaries, but may be incorrect for codimension-3 boundaries 
         *  for graphs with a non-trivial 3 dimensional structure at the time of writing - 9/2/2025
         */
        inline void close();

        /// @brief Number of blocks in the graph
        inline unsigned int numBlocks() const {return m_blocks.size(); }

        /// @brief Number of boundaries block has in direction dir
        /// @param srcBlock a BlockIndex
        /// @param dir a direction specified in the coordinates of srcBlock 
        inline unsigned int numBoundaries(BlockIndex srcBlock, Point& dir) const;
        
        /// @brief Retrieve a vector of arcs associated with boundaries from srcBlock
        /// in the direction dir
        /// @param srcBlock a BlockIndex
        /// @param dir a direction specified in the coordinates of srcBlock
        inline std::vector<MBGraphArc> boundaries(
                BlockIndex srcBlock,
                Point      dir) const;

        /// @brief Returns the list of all directions along which srcBlock is connected to
        /// dstBlock. If two blocks share a codimension-D boundary, this list will contain a single
        /// direction of codimension-D and all of the higher codimensional connections associated with
        /// the shared edges, corners, etc of the boundary between blocks. 
        ///
        /// Intended for internal use only; don't use this unless you know what you are doing.
        ///
        /// @param srcBlock The source block of the implied boundary
        /// @param dstBlock The destination block of the implied boundary
        /// @return The directions associated with the implied boundary. Empty if no boundary exists.
        inline std::vector<Point> fullConnectivity(
                BlockIndex srcBlock,
                BlockIndex dstBlock) const;
        
        /// @brief Returns the lowest codimension direction associated with the boundary between 
        /// two blocks. This direction is guaranteed to be unique since MBGraph does not support
        /// two-block topological cycles. This is the direction associated with the "true" boundary
        /// between blocks.
        ///
        /// @param srcBlock The source block of the implied boundary
        /// @param dstBlock The destination block of the implied boundary
        /// @return The direction (srcBlock coordinates) of the implied boundary. Returns Point::Zeros if no boundary exists.
        /// FIXME: This should return std::optional<Point> with no value if the boundary doesn't exist
        inline Point connectivity(
                BlockIndex srcBlock,
                BlockIndex dstBlock) const;

        /** @brief Returns the direction from dstBlock to srcBlock given the direction from srcBlock to dstBlock.
         *  Fails by assertion if the arc implied by the inputs is not in the graph
         *  FIXME: This should return std::optional<Point> with no value if the arc doesn't exist
         * 
         *  @param srcBlock the source block of the implied boundary
         *  @param dstBlock the destination block of the implied boundary
         *  @param dir the direction associated with an arc from source to destination (srcBlock coordinates)
         *  @return the direction along the same arc as dir from destination to source (dstBlock coordinates)
         */
        inline Point reverseArc(
                BlockIndex srcBlock,
                BlockIndex dstBlock,
                Point dir) const;
        
        /** @brief Given a direction in srcBlock's coordinates, mirror the direction across the boundary with
         *  dstBlock and return the result in dstBlock's coordinates. If R is the coordinate rotation from
         *  srcBlock to dstBlock and dir has normal and tangential components to the boundary n and t respectively,
         *  the output is R(t) - R(n).
         *  
         *  Fails by assertion if there is no boundary between srcBlock and dstBlock
         *  @param srcBlock the source block of the implied boundary
         *  @param dstBlock the destination block of the implied boundary
         *  @param dir a direction in srcBlock's coordinates
         *  @return a direction corresponding to dir being mirrored across the implied boundary
         */
        inline Point mirrorDirAcrossBoundary(
                BlockIndex srcBlock,
                BlockIndex dstBlock,
                Point dir) const;

        /** @brief Attempts to retrieve the coordinate permutation from srcBlock to dstBlock
         *  Fails by assertion if the input blocks don't share a boundary
         * 
         *  @param srcBlock the source block of the implied boundary
         *  @param dstBlock the destination block of the implied boundary
         *  @return the permutation from srcBlock's coordinates to dstBlock's coordinates
         */
         inline CoordPermutation rotation(
            BlockIndex srcBlock,
            BlockIndex dstBlock) const;

        /** @brief Get the block adjacent to srcBlock in the specified codimension-1 direction
         *  Since the direction is codimension-1, the block - should it exist - is guaranteed to be unique
         *  @param srcBlock the block whose adjacency is desired
         *  @param dir the coordinate axis along which to search for an adjacency
         *  @param side the side (high / low) to search for the adjacency on
         */
        inline BlockIndex adjacentBlock(
                BlockIndex      srcBlock,
                unsigned int    dir,
                Side::LoHiSide  side) const;
        
        /** @brief Get the block adjacent to srcBlock in the specified codimension-1 direction
         *  Fails if codim1Dir.codim() != 1 
         *  @param srcBlock a valid block index
         *  @param codim1Dir a codimension-1 direction
        */
        inline BlockIndex adjacentBlock(
                BlockIndex    srcBlock,
                Point         codim1Dir) const;
        
        /** @brief Get a vector of all blocks adjacent to this block in direction dir
         *  
         *  @param srcBlock a valid block index
         *  @param dir a direction
         */
        inline std::set<BlockIndex> adjacentBlocks(
                BlockIndex srcBlock,
                Point dir) const;

        /** @brief Equality operator. Equality is determined by value */
        inline bool operator==(const MBGraph& rhs) const;

        /** @brief Determine if the block - direction pair represents the boundary of a block
         *  This function returns true if dir points to a "triple point" region even though the
         *  graph has no arcs from the node associated with srcBlock in the direction dir
         *  
         *  @param srcBlock a valid block index
         *  @param dir a direction
         */
        inline bool isBlockBoundary(BlockIndex srcBlock, Point dir) const;
      
        
        /** @brief returns true if the specified block boundary corresponds to a "triple point".
         *  In this context, a "triple point" is a codimension-2 block boundary where exactly
         *  three blocks meet. This function is used to infer corrections that must be made 
         *  when building interpolation operators in such regions.
         * 
         *  Example: A "Triple Point" Configuration in 2D
         *        \
         *         \ Block 2
         *  Block 0 >-------
         *         / Block 1
         *        /
         */
        inline bool isTriplePoint(
                BlockIndex  srcBlock,
                Point        dir) const;

        /** @brief return a list of codimension-1 direction pairs which represent "triple point circuits"
         *  two codimension-1 directions represent a circuit relative to srcBlock if srcBlock has a
         *  codimension-1 boundary in those directions with two other blocks which in turn have a codimenison-1
         *  boundary which lies in the same logical "plane" as the two directions defining the circuit.
         *  
         *  This function is not recommended for public use.
         * 
         *  If DIM==2, there will be at most one circuit for a given block. In higher dimensions, there will in
         *  general be many, and will result in the function being called recursively.
         * 
         *  At the time of writing, there are no dependencies in Proto which use this function, but it is likely
         *  that it will be useful for more complex graphs. - Chris Bozhart - 9/2/2025
         * 
         *  @param srcBlock the reference block from which the search for triple point circuits is executed
         *  @param dir the direction in which triple point circuits are searched
        */
        std::vector<std::pair<Point, Point>>
        getTriplePointCircuits(
            BlockIndex srcBlock,
            Point dir) const;

        /** @brief returns true if there are no blocks adjacent to srcBlock in the direction dir. This situation
         *  is assumed to correspond to a "domain boundary"
         *  @param srcBlock the reference block from which domain boundaries are searched for
         *  @param dir the direction which constitutes the query for a domain boundary
         */
        inline bool isDomainBoundary(
                BlockIndex  srcBlock,
                Point       dir) const;

        /// @brief prints all of the MBGraphNode objects in this graph (one for each block)
        /** which in turn prints all of the arcs extending from each block
        */
        inline void print() const;
        

        
    private:
        
        /** @brief internal function called inside of "close" adjusts errors in the implicitly constructed
         *  coordinate permutations between blocks with high codimensional connectivity. 
         */
        inline void fixRotations();

        /** @brief obtain a reference to the coordinate permuation object stored in the graph. Used
         *  by "close" to fix internal boundary rotations, and not intended for public use. If you think you 
         *  need to use this externally, there is probably a bug in MBGraph::close that needs addressing instead
         */
        inline CoordPermutation& rotation(
                BlockIndex srcBlock,
                Point dir,
                BlockIndex dstBlock) const;
    
        /** @brief obtain a reference to the coordinate permuation object stored in the graph. Used
         *  by "close" to fix internal boundary rotations, and not intended for public use. If you think you 
         *  need to use this externally, there is probably a bug in MBGraph::close that needs addressing instead
         */
        inline CoordPermutation& rotation(
                BlockIndex block,
                Point dir) const;

        /** @brief for a given block - direction pair, return all of the codimension-2 boundaries 
         *  that are projections of dir which also represent triple point boundaries.
         * 
         *  This function is highly experimental and not well tested.
         *  This function only has true utility for DIM >= 3. For an illustrative example of when
         *  This function has value, consider a graph of 4 blocks connected in such a way that
         *  the arcs for a tetrahedron. It is completely unused at the time of writing, but
         *  may come in handy in the future.
         *  - Chris Bozhart - 8/28/2025
         */
        inline std::vector<Point> codim2TriplePointBoundaries(
                BlockIndex srcBlock,
                Point dir) const;


        bool m_closed; ///< Has "close" been called on this
        std::vector<std::shared_ptr<MBGraphNode>> m_blocks; ///< The nodes of this graph
    };
#include "implem/Proto_MBGraphImplem.H"
} // end namespace Proto
#endif // end include guard
