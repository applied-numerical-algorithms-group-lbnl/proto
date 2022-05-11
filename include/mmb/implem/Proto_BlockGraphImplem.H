
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
    m_index = a_index;
    m_boundaries.clear();
    m_rotations.clear();
}

bool BlockGraphNode::addBoundary(
    Point                   a_dir,
    std::shared_ptr<BlockGraphNode>  a_block,
    CoordPermutation&       a_rotation)
{
    PROTO_ASSERT(codim(a_dir) == 1,
        "BlockGraphNode::addBoundary | Error: \
        Permutations may only be specified for boundaries of codimension 1.");
    PROTO_ASSERT(m_boundaries.find(a_dir) == m_boundaries.end(),
        "BlockGraphNode::addBoundary | Error: \
        Attempting to add a duplicate codimension 1 boundary to BlockGraphNode %u.",
        m_index);
    
    m_rotations[a_dir] = a_rotation;
    return addBoundary(a_dir, a_block);
}

bool BlockGraphNode::addBoundary(
        Point                   a_dir,
        std::shared_ptr<BlockGraphNode>  a_block)
{
    if (a_block->id() != this->id()) {return false;}
    m_boundaries[a_dir].insert(a_block);
    return true;
}

const std::set<std::shared_ptr<BlockGraphNode>>& BlockGraphNode::boundaries(Point& a_dir) const
{
    PROTO_ASSERT(numBoundaries(a_dir) > 0,
        "BlockGraphNode::boundaries | Error: Tried to access boundaries that don't exist.");
    return (const_cast<BlockGraphNode*>(this))->m_boundaries[a_dir];
}

unsigned int BlockGraphNode::numBoundaries(Point& a_dir) const
{
    if (m_boundaries.find(a_dir) == m_boundaries.end())
    {
        return 0;
    } else {
        return (const_cast<BlockGraphNode*>(this))->m_boundaries[a_dir].size();
    }
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
    
    srcBlock->addBoundary(dir_01, dstBlock, a_rotation);
    auto inverseRotation = a_rotation.inverse();
    dstBlock->addBoundary(dir_10, srcBlock, inverseRotation);

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
    
        if (srcBlock->numBoundaries(dir_i_01) > 0)
        {
            for (auto block : srcBlock->boundaries(dir_i_01))
            {
                dstBlock->addBoundary(dir_i_10, block);
            }
        }
        if (dstBlock->numBoundaries(dir_i_10) > 0)
        {
            for (auto block : dstBlock->boundaries(dir_i_10))
            {
                srcBlock->addBoundary(dir_i_01, block);
            }
        }
        srcBlock->addBoundary(dir_i_01, dstBlock);
        dstBlock->addBoundary(dir_i_10, srcBlock);
    }
}

void BlockGraph::closeCircuit(
        unsigned int        a_srcBlock,
        unsigned int        a_dstBlock)
{
#if DIM >= 2
    auto srcBlock = m_blocks[a_srcBlock];
    auto dstBlock = m_blocks[a_dstBlock];
    // there must be exactly 1 codim-2 boundary between src and dst
#endif
}

const BlockGraphNode& BlockGraph::operator[](unsigned int a_index) const
{
    return *m_blocks[a_index];
}

unsigned int BlockGraph::numBoundaries(unsigned int a_index, Point& a_dir) const
{
    PROTO_ASSERT(a_index < size(),
            "BlockGraph::numBoundaries | Error: Block index %u is out of bounds (size = %u).",
            a_index, size());
    return (*this)[a_index].numBoundaries(a_dir); 
}

std::set<Point> BlockGraph::connectivity(unsigned int a_srcBlock, unsigned int a_dstBlock) const
{
    auto& srcBlock = *m_blocks[a_srcBlock];
    std::set<Point> ret;
    if (a_srcBlock == a_dstBlock) { return ret; }
    Box K = Box::Kernel(1);
    for (auto dir : K)
    {
        if (numBoundaries(a_srcBlock, dir) > 0)
        {
            auto bounds = srcBlock.boundaries(dir);
            for (auto block : bounds)
            {
                if (block->id() == a_dstBlock) { ret.insert(dir); }
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
    auto bounds = m_blocks[a_srcBlock]->boundaries(dir);
    PROTO_ASSERT(bounds.size() == 1,
            "BlockGraph::adjacent | Error: \
            Data corruption. Block %u has multiple blocks bounding the same face.",
            a_srcBlock);
    return (*bounds.begin())->id();
}

void BlockGraph::print() const
{
    for (auto block : m_blocks)
    {
        pout() << "Block " << block->id() << std::endl;
        Box K = Box::Kernel(1);
        for (auto dir : K)
        {
            if (block->numBoundaries(dir) > 0)
            {
                pout() << "\tDir: " << dir << " (codim = " << codim(dir) << ")" << std::endl;
                for (auto bound : block->boundaries(dir))
                {
                    pout() << "\t\tBound: Block " << bound->id() << std::endl;
                }
            }
        }
    }
}
