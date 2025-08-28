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

        inline void print(std::ostream& stream = std::cout) { stream << str() << std::endl; }

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



        /// Get Adjacent Block
        /** Obtain the block adjacent to a specified block in a specified codimension-1
         *  direction */
        inline BlockIndex adjacent(
                BlockIndex    srcBlock,
                unsigned int    dir,
                Side::LoHiSide  side) const;
        
        /// Get Adjacent Block (Overload)
        /** This function is the same as the other version of adjacent. This function will fail
         *  if a codimension-2+ point is supplied. */
        inline BlockIndex adjacent(
                BlockIndex    srcBlock,
                Point           dir) const;
        
        inline std::set<BlockIndex> adjacentBlocks(
                BlockIndex srcBlock,
                Point dir) const;

        /// Equality Operator
        inline bool operator==(const MBGraph& rhs) const;

        /// Detect Block Boundary
        /** Returns true if the specified inputs correspond to a boundary between two blocks. */
        inline bool isBlockBoundary(BlockIndex srcBlock, Point dir) const;
      
        /// Return list of codimension-2 triple points
        inline std::vector<Point> triplePoints(
                BlockIndex srcBlock,
                Point dir) const;

        /// Detect Triple Point
        /** Returns true if the specified block boundary corresponds to a "triple point".
         *  In this context, a "triple point" is a codimension-2 block boundary where exactly
         *  three blocks meet. This function is used to infer corrections that must be made 
         *  when building interpolation operators in such regions. */
        inline bool isTriplePoint(
                BlockIndex  srcBlock,
                Point        dir) const;

        std::vector<std::pair<Point, Point>>
        getTriplePointCircuits(
            BlockIndex srcBlock,
            Point dir) const;

        inline bool isDomainBoundary(
                BlockIndex  srcBlock,
                Point       dir) const;

        /// Print
        inline void print() const;
        
        /// Close Graph
        /** Executes postprocessing steps once all boundaries have been defined. This function
         *  Also checks for various topological errors in the graph construction. */
        inline void close();
        
    private:
    
        inline void fixRotations();

        inline CoordPermutation& rotation(
                BlockIndex srcBlock,
                Point dir,
                BlockIndex dstBlock) const;

        inline CoordPermutation& rotation(
                BlockIndex block,
                Point dir) const;

        bool m_closed;
        std::vector<std::shared_ptr<MBGraphNode>> m_blocks;
    };
#include "implem/Proto_MBGraphImplem.H"
} // end namespace Proto
#endif // end include guard
