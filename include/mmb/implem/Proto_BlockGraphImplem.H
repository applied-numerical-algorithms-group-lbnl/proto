
unsigned int codim(Point a_dir)
{
    int numZeros = 0;
    for (int d = 0; d < DIM; d++)
    {
        if (a_dir[d] != 0) { numZeros++; }
    }
    return numZeros;
}

std::vector<Point> codimDirs(unsigned int a_codim)
{
    PROTO_ASSERT(a_codim <= DIM,
            "codimDirs | Error: Invalid codimension value: %u", a_codim);
    std::vector<Point> dirs;
    Box K = Box::Kernel(1);
    for (auto dir : K)
    {
        if (codim(dir) == a_codim)
        {
            dirs.push_back(dir);
        }
    }
    return dirs;
}

//=================================================================================================
// BLOCK GRAPH NODE

BlockGraphNode::BlockGraphNode(unsigned int a_index)
{
    index = a_index;
    boundaries.clear();
    Box K = Box::Kernel(1);
    for (auto dir : K)
    {
        if (dir == Point::Zeros()) { continue; }
        boundaries[dir].clear();
    }
    rotations.clear();
}
void BlockGraphNode::addArc(Point a_dir, ArcPair& a_arc)
{
    // ignore self arcs
    if (a_arc.first->index == this->index) { return; }
    for (auto arc : boundaries[a_dir])
    {
        // ignore duplicate arcs
        if (arc == a_arc) { return; }
    }
    boundaries[a_dir].push_back(a_arc);
}

//=================================================================================================
// BLOCK GRAPH

BlockGraph::BlockGraph(unsigned int a_numBlocks)
{
    m_blocks.resize(a_numBlocks);
    for (int ii = 0; ii < a_numBlocks; ii++)
    {
        m_blocks[ii] = std::make_shared<BlockGraphNode>(ii);
    }
}

void BlockGraph::addBoundary(
        unsigned int        a_srcBlock,
        unsigned int        a_dstBlock,
        unsigned int        a_dir,
        Side::LoHiSide      a_side,
        CoordPermutation&   a_rotation)
{
    // TODO: Remove this restriction if we allow for periodic boundaries
    PROTO_ASSERT(a_srcBlock != a_dstBlock,
        "BlockGraph::addBoundary | Error: src and dst block indices must be different.");
    PROTO_ASSERT(a_srcBlock < size(),
        "BlockGraph::addBoundary | Error: srcBlock index %u is out of bounds (size = %u).",
        a_srcBlock, size());
    PROTO_ASSERT(a_dstBlock < size(),
        "BlockGraph::addBoundary | Error: dstBlock index %u is out of bounds (size = %u).",
        a_dstBlock, size());

    auto srcBlock = m_blocks[a_srcBlock];
    auto dstBlock = m_blocks[a_dstBlock];
   
    Point dir_01 = Point::Basis(a_dir, a_side);
    Point dir_10 = a_rotation(-dir_01);
    
    PROTO_ASSERT(srcBlock->boundaries[dir_01].size() == 0,
            "BlockGraph::addBoundary | Error: Adding redundant codimension-1 boundary to block %u",
            a_srcBlock);
    PROTO_ASSERT(dstBlock->boundaries[dir_10].size() == 0,
            "BlockGraph::addBoundary | Error: Adding redundant codimension-1 boundary to block %u",
            a_dstBlock);

    auto inverseRotation = a_rotation.inverse();
    
    srcBlock->boundaries[dir_01].push_back(ArcPair(dstBlock, dir_10));
    srcBlock->rotations[dir_01] = a_rotation;
    dstBlock->boundaries[dir_10].push_back(ArcPair(srcBlock, dir_01));
    dstBlock->rotations[dir_10] = inverseRotation;

    // Add all boundaries of codimension > 1
    Box dirBox_01 = Box::Kernel(1).grow(a_dir, -1);
    for (auto biter = dirBox_01.begin(); biter.ok(); ++biter)
    {
        if (*biter == Point::Zeros()) { continue; }
        Point dir_i_01 = *biter + dir_01;
        Point dir_i_10 = a_rotation(*biter) + dir_10;
        
        PROTO_ASSERT(codim(dir_i_01) > 1,
            "BlockGraph::addBoundary | boundary has invalid codimension");
        PROTO_ASSERT(codim(dir_i_10) > 1,
            "BlockGraph::addBoundary | boundary has invalid codimension");
    
        for (auto arc : srcBlock->boundaries[dir_i_01])
        {
            ArcPair toDst(dstBlock, arc.second);
            arc.first->addArc(arc.second, toDst);
            ArcPair fromDst(arc.first, dir_i_10);
            dstBlock->addArc(dir_i_10, fromDst);
        }
        for (auto arc : dstBlock->boundaries[dir_i_10])
        {
            ArcPair toSrc(srcBlock, arc.second);
            arc.first->addArc(arc.second, toSrc);
            ArcPair fromSrc(arc.first, dir_i_01);
            srcBlock->addArc(dir_i_01, fromSrc);
        }
        ArcPair srcToDst(dstBlock, dir_i_10);
        srcBlock->addArc(dir_i_01, srcToDst);
        ArcPair dstToSrc(srcBlock, dir_i_01);
        dstBlock->addArc(dir_i_10, dstToSrc);
    }
}

void BlockGraph::closeCircuit(
        unsigned int        a_srcBlock,
        unsigned int        a_dstBlock)
{
#if DIM >= 2
    auto connections = connectivity(a_srcBlock, a_dstBlock);
    std::vector<Point> codim2Bounds;
    for (auto conn : connections)
    {
        PROTO_ASSERT(codim(conn) > 1,
            "BlockGraph::closeCircuit | Error: \
            Blocks %u and %u already share a codimension 1 boundary. There is no circuit \
            to close.", a_srcBlock, a_dstBlock);

        if (codim(conn) == 2) { codim2Bounds.push_back(conn); }
    }
    PROTO_ASSERT(codim2Bounds.size() == 1,
            "BlockGraph::closeCircuit | Error: \
            Could not close circuit between blocks %u and %u. \
            A circuit can only be closed between two blocks if they share exactly \
            one boundary of codimensioon 2.", a_srcBlock, a_dstBlock);

    Point dir_01 = codim2Bounds[0];
#endif
}

unsigned int BlockGraph::numBoundaries(unsigned int a_index, Point& a_dir) const
{
    PROTO_ASSERT(a_index < size(),
            "BlockGraph::numBoundaries | Error: Block index %u is out of bounds (size = %u).",
            a_index, size());
    return m_blocks[a_index]->boundaries[a_dir].size();
}

std::set<Point> BlockGraph::connectivity(unsigned int a_srcBlock, unsigned int a_dstBlock) const
{
    auto& srcBlock = *m_blocks[a_srcBlock];
    std::set<Point> ret;
    Box K = Box::Kernel(1);
    for (auto dir : K)
    {
        if (dir == Point::Zeros()) { continue; }
        for (auto arc : srcBlock.boundaries[dir])
        {
            if (arc.first->index == a_dstBlock)
            {
                ret.insert(dir);
            }
        }
    }
    return ret; 
}

unsigned int BlockGraph::adjacent(
        unsigned int    a_srcBlock,
        unsigned int    a_dir,
        Side::LoHiSide  a_side) const
{
    Point dir = Point::Basis(a_dir, a_side);
    if (numBoundaries(a_srcBlock, dir) == 0) { return size(); }
    auto arcs = m_blocks[a_srcBlock]->boundaries[dir];
    PROTO_ASSERT(arcs.size() == 1,
            "BlockGraph::adjacent | Error: \
            Data corruption. Block %u has multiple blocks bounding the same face.",
            a_srcBlock);
    return arcs[0].first->index;
}

void BlockGraph::print() const
{
    for (auto block : m_blocks)
    {
        pout() << "Block " << block->index << std::endl;
        Box K = Box::Kernel(1);
        for (auto dir : K)
        {
            pout() << "\tDir: " << dir << " (codim = " << codim(dir) << ")" << std::endl;
            for (auto arc : block->boundaries[dir])
            {
                pout() << "\t\tBound: Block " << arc.first->index << std::endl;
            }
        }
    }
}
