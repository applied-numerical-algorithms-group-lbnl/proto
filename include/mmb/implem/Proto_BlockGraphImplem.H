
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

std::vector<Point> getBasisVectors(const Point& a_pt)
{
    std::vector<Point> basis;
    for (int dir = 0; dir < DIM; dir++)
    {
        if (a_pt[dir] != 0)
        {
            basis.push_back(Point::Basis(dir)*a_pt[dir]);
        }
    }
    PROTO_ASSERT(codim(a_pt) == basis.size(),
            "getBasisVector | Error: Data corruption.");
    return basis;
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
    auto dst = a_arc.first->index;
    if (dst == this->index) { return; }
    for (auto arc : boundaries[a_dir])
    {
        // ignore duplicate arcs
        if (arc.first->index == dst) { return; }
    }
    //std::cout << "Adding Arc: " << a_dir << " --> " << a_arc.first->index << " (<-- " << a_arc.second << ")" << std::endl;
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
    addBoundary(a_srcBlock, a_dstBlock, Point::Basis(a_dir, a_side), a_rotation);
}
void BlockGraph::addBoundary(
        unsigned int        a_srcBlock,
        unsigned int        a_dstBlock,
        Point               a_dir,
        CoordPermutation&   a_rotation)
{
    
    //pout() << "Adding Boundary | src: " << a_srcBlock << " | dst: " << a_dstBlock << " | dir: " << a_dir << std::endl;
    //std::cout << "Adding Boundary | src: " << a_srcBlock << " | dst: " << a_dstBlock << " | dir: " << a_dir << std::endl;
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
   
    Point dir_01 = a_dir;
    Point dir_10 = a_rotation(-dir_01);
   
    //if (srcBlock->boundaries[dir_01].size() != 0) { return; } //redundant boundary
    //if (srcBlock->boundaries[dir_10].size() != 0) { return; } //redundant boundary
    PROTO_ASSERT(srcBlock->boundaries[dir_01].size() == 0,
            "BlockGraph::addBoundary | Error: Adding redundant codimension-1 boundary to block %u",
            a_srcBlock);
    PROTO_ASSERT(dstBlock->boundaries[dir_10].size() == 0,
            "BlockGraph::addBoundary | Error: Adding redundant codimension-1 boundary to block %u",
            a_dstBlock);
    auto inverseRotation = a_rotation.inverse();
   
    ArcPair arc_01(dstBlock, dir_10);
    srcBlock->addArc(dir_01, arc_01);
    ArcPair arc_10(srcBlock, dir_01);
    dstBlock->addArc(dir_10, arc_10);

    srcBlock->rotations[dir_01] = a_rotation;
    dstBlock->rotations[dir_10] = inverseRotation;

    // Add all boundaries of codimension > 1
    unsigned int normDir;
    for (int d = 0; d < DIM; d++) { if (a_dir[d] != 0) { normDir = d; break; } }
    Box dirBox_01 = Box::Kernel(1).grow(normDir, -1);
    for (auto biter = dirBox_01.begin(); biter.ok(); ++biter)
    {
        if (*biter == Point::Zeros()) { continue; }
        Point dir_i_01 = *biter + dir_01;
        Point dir_i_10 = a_rotation(*biter) + dir_10;
        
        PROTO_ASSERT(codim(dir_i_01) > 1,
            "BlockGraph::addBoundary | boundary has invalid codimension");
        PROTO_ASSERT(codim(dir_i_10) > 1,
            "BlockGraph::addBoundary | boundary has invalid codimension");
        
        ArcPair toDst(dstBlock, dir_i_10);
        ArcPair toSrc(srcBlock, dir_i_01);
        for (auto arc : srcBlock->boundaries[dir_i_01])
        {
            /*
            ArcPair toDst(dstBlock, arc.second);
            arc.first->addArc(arc.second, toDst);
            ArcPair fromDst(arc.first, dir_i_10);
            dstBlock->addArc(dir_i_10, fromDst);
            */
            arc.first->addArc(arc.second, toDst);
            dstBlock->addArc(dir_i_10, arc);
        }
        for (auto arc : dstBlock->boundaries[dir_i_10])
        {
            /*
            ArcPair toSrc(srcBlock, arc.second);
            arc.first->addArc(arc.second, toSrc);
            ArcPair fromSrc(arc.first, dir_i_01);
            srcBlock->addArc(dir_i_01, fromSrc);
            */
            arc.first->addArc(arc.second, toSrc);
            srcBlock->addArc(dir_i_01, arc);
        }
        ArcPair srcToDst(dstBlock, dir_i_10);
        srcBlock->addArc(dir_i_01, srcToDst);
        ArcPair dstToSrc(srcBlock, dir_i_01);
        dstBlock->addArc(dir_i_10, dstToSrc);
    }
    //print();
    //pout() << "====================================================================================" << std::endl;
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

    auto srcBlock = m_blocks[a_srcBlock];
    auto dstBlock = m_blocks[a_dstBlock];

    auto edgeArcs = srcBlock->boundaries[codim2Bounds[0]];
    Point dir_01;
    for (auto dir : getBasisVectors(codim2Bounds[0]))
    {
        if (adjacent(a_srcBlock, dir) == size())
        {
            dir_01 = dir;
        }
    }

    CoordPermutation circuitRotation;
    unsigned int currBlock = a_srcBlock;
    unsigned int prevBlock = size(); //null value
    int iter = 0;
    while (currBlock != a_dstBlock)
    {
        pout() << "Current Block: " << currBlock << std::endl;
        pout() << "Current permutation: " << std::endl;
        circuitRotation.print();
        // get the dir corresponding to the codim-2 boundary in currBlock coords
        Point edgeDir;
        if (currBlock == a_srcBlock) { edgeDir = codim2Bounds[0]; }
        else {
            for (auto arc : edgeArcs)
            {
                if (arc.first->index == currBlock)
                {
                    edgeDir = arc.second;
                    break;
                }
            }
        }
        // get the codim-1 boundaries which intersect at edgeDir
        auto adjDirs = getBasisVectors(edgeDir);
        // advance around the circuit
        for (auto adjDir : adjDirs)
        {
            unsigned int adjBlock = adjacent(currBlock, adjDir);
            if (adjBlock != prevBlock)
            {
                auto adjRot = m_blocks[currBlock]->rotations[adjDir];
                pout() << "adjactent rotation: " << std::endl;
                adjRot.print();
                circuitRotation = circuitRotation*adjRot;
                prevBlock = currBlock;
                currBlock = adjBlock;
            }
        }
        if (currBlock == a_srcBlock)
        {
            // circuit was already closed
            return;
        }
        iter++;
        PROTO_ASSERT(iter < size(),
            "BlockGraph::closeCircuit | Error: \
            Detected a circuit which is larger than the size of the graph.");
    }
    CoordPermutation testRot = circuitRotation * (srcBlock->rotations[dir_01].inverse());
    PROTO_ASSERT(testRot == CoordPermutation::identity(),
            "BlockGraph::closeCircuit | Error: \
            CoordPermutation around the circuit is not identity.");
    std::cout << "Closing circuit" << std::endl;
    addBoundary(a_srcBlock, a_dstBlock, dir_01, circuitRotation);
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
    return adjacent(a_srcBlock, Point::Basis(a_dir, a_side));
}
unsigned int BlockGraph::adjacent(
        unsigned int    a_srcBlock,
        Point a_dir) const
{
    if (numBoundaries(a_srcBlock, a_dir) == 0) { return size(); }
    auto arcs = m_blocks[a_srcBlock]->boundaries[a_dir];
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
                pout() << "\t\tBound: Block " << arc.first->index << " | Local Dir: " << arc.second << std::endl;
            }
        }
    }
}
