unsigned int codim(Point a_dir)
{
    int numZeros = 0;
    for (int d = 0; d < DIM; d++)
    {
        if (a_dir[d] != 0) { numZeros++; }
    }
    return numZeros;
}

//=================================================================================================
//  BLOCK


Block::Block(unsigned int a_index)
{
    m_index = a_index;
    m_boundaries.clear();
    m_rotations.clear();
}

bool Block::addBoundary(
    Point                   a_dir,
    std::shared_ptr<Block>  a_block,
    CoordPermutation        a_rotation)
{
    PROTO_ASSERT(codim(a_dir) == 1,
        "Block::addBoundary | Error: \
        Permutations may only be specified for boundaries of codimension 1.");
    PROTO_ASSERT(m_boundaries.find(a_dir) == m_boundaries.end(),
        "Block::addBoundary | Error: \
        Attempting to add a duplicate codimension 1 boundary to Block %u.",
        m_index);
    PROTO_ASSERT(a_block->id() != this->id(),
        "Block::addBoundary | Error: Blocks are not allowed to have boundaries to themselves.");
    
    m_rotations[a_dir] = a_rotation;
    addBoundary(a_dir, a_block);
}

bool Block::addBoundary(
        Point                   a_dir,
        std::shared_ptr<Block>  a_block)
{
    PROTO_ASSERT(m_boundaries.find(a_dir) == m_boundaries.end(),
        "Block::addBoundary | Error: \
        Attempting to add a duplicate codimension 1 boundary to Block %u.",
        m_index);
    PROTO_ASSERT(a_block->id() != this->id(),
        "Block::addBoundary | Error: Blocks are not allowed to have boundaries to themselves.");
    m_boundaries[a_dir] = a_block;
}

const std::vector<std::shared_ptr<Block>>& Block::operator[](Point& a_dir) const
{
    PROTO_ASSERT(numBoundaries(a_dir) > 0,
        "Block::operator[] | Error: Tried to access boundaries that don't exist.");
    return m_boundaries[a_dir];
}

unsigned int Block::numBoundaries(Point& a_dir) const
{
    if (m_boundaries.find(a_dir) == m_boundaries.end())
    {
        return 0;
    } else {
        return m_boundaries[a_dir].size();
    }
}
//=================================================================================================
//  BLOCK GRAPH

BlockGraph::BlockGraph(unsigned int a_numBlocks, Point a_blockSize)
{
    m_blockSize = a_blockSize;
    m_blocks.resize(a_numBlocks);
    for (int ii = 0; ii < a_numBlocks; ii++)
    {
        m_blocks[ii] = std::make_shared<Block>(ii, a_blockSize);
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
    
    Block& srcBlock = *m_blocks[a_srcBlock];
    Block& dstBlock = *m_blocks[a_dstBlock];

    Point dir_01(a_dir, a_side);
    Point dir_10 = a_rotation(-dir_01);

    srcBlock.addBoundary(dir_01, dstBlock, a_rotation);
    dstBlock.addBoundary(dir_10, srcBlock, a_rotation.inverse());

    // Add all boundaries of codimension > 1
    Box dirBox_01 = Box::Kernel(3).grow(a_dir, -1);
    for (auto biter = dirBox_01.begin(); biter.ok(); ++biter)
    {
        if (*biter == Point::Zeros()) { continue; }
        Point dir_i_01 = *biter + dir_01;
        Point dir_i_10 = a_rotation(*biter) + dir_10;

        PROTO_ASSERT(codim(dir_i_01) > 1,
            "BlockGraph::addBoundary | boundary has invalid codimension");
        PROTO_ASSERT(codim(dir_i_10) > 1,
            "BlockGraph::addBoundary | boundary has invalid codimension");

        srcBlock.addBoundary(dir_i_01, dstBlock);
        dstBlock.addBoundary(dir_i_10, srcBlock);
    }
}


