std::vector<Point> codimDirs(unsigned int a_codim)
{
    PROTO_ASSERT(a_codim <= DIM,
            "codimDirs | Error: Invalid codimension value: %u", a_codim);
    std::vector<Point> dirs;
    Box K = Box::Kernel(1);
    for (auto dir : K)
    {
        if (dir.codim() == a_codim)
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
    PROTO_ASSERT(a_pt.codim() == basis.size(),
            "getBasisVector | Error: Data corruption.");
    return basis;
}

//=================================================================================================
// BLOCK GRAPH NODE

MBGraphNode::MBGraphNode(unsigned int a_index)
{
    index = a_index;
    arcs.clear();
}
void MBGraphNode::addArc(unsigned int a_block,
    Point a_toDst, Point a_fromDst)
{
    if (a_block == index) { return; }
    MBGraphArc arc(index, a_block, a_toDst, a_fromDst); 
    for (auto arc : arcs[a_toDst])
    {
        if (arc.dstBlock == a_block) {return; } //redundant arc
    }
    arcs[a_toDst].push_back(arc);
}

CoordPermutation& MBGraphNode::rotation(Point a_dir, unsigned int a_block)
{
    for (auto& bound : arcs[a_dir])
    {
        if (bound.dstBlock == a_block)
        {
            return bound.R;
        }
    }
    MayDay<void>::Error("MBGraphNode::rotation | Error: Data corruption.");
}

//=================================================================================================
// BLOCK GRAPH

MBGraph::MBGraph(unsigned int a_numBlocks)
{
    m_closed = false;
    m_blocks.resize(a_numBlocks);
    for (int ii = 0; ii < a_numBlocks; ii++)
    {
        m_blocks[ii] = std::make_shared<MBGraphNode>(ii);
    }
}

void MBGraph::defineBoundary(
        unsigned int        a_srcBlock,
        unsigned int        a_dstBlock,
        unsigned int        a_dir,
        Side::LoHiSide      a_side,
        CoordPermutation&   a_rotation)
{
    defineBoundary(a_srcBlock, a_dstBlock, Point::Basis(a_dir, a_side), a_rotation);
}

void MBGraph::fixRotations()
{
    PROTO_ASSERT(!m_closed,
        "MBGraph::fixRotations | Error: Graph definition has been completed and is locked");
    for (unsigned int b0 = 0; b0 < size(); b0++)
    for (auto d0 : codimDirs(2))
    {
        auto dirs = d0.split();
        Point f0 = dirs[0];
        Point g0 = dirs[1];
       
        bool done = false;
        bool reverse = false;
        while (!done)
        {
            Point d = d0;
            Point f = f0;
            Point g = g0;
            unsigned int b1 = b0;
            auto b2 = adjacent(b0,f0);
            while (b2 < size() && b2 != b0)
            {
                d = reverseArc(b1,b2,d);
                g = reverseArc(b1,b2,f);
                f = d - g;
                if (connectivity(b0,b2).codim() == 2)
                {
                    CoordPermutation R02{{f0,-f},{g0,-g}};
                    
                    m_blocks[b0]->rotation(d0,b2) = R02;
                }

                b1 = b2;
                b2 = adjacent(b2,f);
            }
            if (b2 == b0 || reverse)
            {
                done = true;
            } else {
                reverse = true;
            }
        }
    }
}

void MBGraph::fixRotationsOld()
{
    pout() << "MBGraph::fixRotations" << std::endl;
    for (int c = 1; c < DIM; c++)
    for (unsigned int b0 = 0; b0 < size(); b0++)
    for (auto d01 : codimDirs(c))
    for (auto bound_01 : boundaries(b0, d01))
    {
        unsigned int b1 = bound_01.dstBlock;
        for (auto d02 : codimDirs(c+1))
        {
            if (d01.dot(d02) <= 0) { continue; } //d02 is not in d01
            auto bounds_02 = boundaries(b0, d02);
            if (bounds_02.size() == 0) { continue; } //there are no blocks in d02
            auto R01 = rotation(b0, d01, b1);
            auto d12 = R01(d02 - d01);
            PROTO_ASSERT(d12.codim() == 1, "MBGraph::fixRotations | Error: Data corruption.");
            unsigned int b2 = adjacent(b1,d12);
            if (b2 >= size()) { continue; } // b1 is not part of a valid circuit
            auto R12 = rotation(b1,d12,b2);
            auto R02 = R01*R12;
            m_blocks[b0]->rotation(d02,b2) = R02;
        }
    }
}

void MBGraph::defineBoundary(
        unsigned int        a_srcBlock,
        unsigned int        a_dstBlock,
        Point               a_dir,
        CoordPermutation&   a_rotation)
{
    PROTO_ASSERT(!m_closed,
        "MBGraph::defineBoundary | Error: Graph definition has been completed and is locked");
    // TODO: Remove this restriction if we allow for periodic boundaries
    PROTO_ASSERT(a_srcBlock != a_dstBlock,
        "MBGraph::defineBoundary | Error: src and dst block indices must be different.");
    PROTO_ASSERT(a_srcBlock < size(),
        "MBGraph::defineBoundary | Error: srcBlock index %u is out of bounds (size = %u).",
        a_srcBlock, size());
    PROTO_ASSERT(a_dstBlock < size(),
        "MBGraph::defineBoundary | Error: dstBlock index %u is out of bounds (size = %u).",
        a_dstBlock, size());
    PROTO_ASSERT(a_dir.codim() == 1,
        "MBGraph::defineBoundary | Error: Only codimension-1 boundaries may be manually defined.");

    auto srcNode = m_blocks[a_srcBlock];
    auto dstNode = m_blocks[a_dstBlock];
   
    Point dir_01 = a_dir;
    Point dir_10 = a_rotation(-dir_01);
   
    PROTO_ASSERT(srcNode->arcs[dir_01].size() == 0,
            "MBGraph::defineBoundary | Error: Adding redundant codimension-1 boundary to block %u",
            a_srcBlock);
    PROTO_ASSERT(dstNode->arcs[dir_10].size() == 0,
            "MBGraph::defineBoundary | Error: Adding redundant codimension-1 boundary to block %u",
            a_dstBlock);
    auto inverseRotation = a_rotation.inverse();
   
    srcNode->addArc(a_dstBlock, dir_01, dir_10);
    dstNode->addArc(a_srcBlock, dir_10, dir_01);
    srcNode->rotation(dir_01, a_dstBlock) = a_rotation;
    dstNode->rotation(dir_10, a_srcBlock) = inverseRotation;
    // Add all boundaries of codimension > 1
    Point normDir = a_dir*a_dir;
    Box dirBox_01 = Box::Kernel(1).grow(-normDir);
    for (auto biter : dirBox_01)
    {
        if (biter == Point::Zeros()) { continue; }
        Point dir_i_01 = biter + dir_01;
        Point dir_i_10 = a_rotation(biter) + dir_10;
        
        PROTO_ASSERT(dir_i_01.codim() > 1,
            "MBGraph::defineBoundary | boundary has invalid codimension");
        PROTO_ASSERT(dir_i_10.codim() > 1,
            "MBGraph::defineBoundary | boundary has invalid codimension");
        
        for (auto arc : srcNode->arcs[dir_i_01])
        {
            auto otherBlock = arc.dstBlock;
            auto otherNode = m_blocks[otherBlock];
            otherNode->addArc(a_dstBlock, arc.dstToSrc, dir_i_10);
            dstNode->addArc(otherBlock, dir_i_10, arc.dstToSrc);
        }
        for (auto arc : dstNode->arcs[dir_i_10])
        {
            auto otherBlock = arc.dstBlock;
            auto otherNode = m_blocks[otherBlock];
            otherNode->addArc(a_srcBlock, arc.dstToSrc, dir_i_01);
            srcNode->addArc(otherBlock, dir_i_01, arc.dstToSrc);
        }
        srcNode->addArc(a_dstBlock, dir_i_01, dir_i_10);
        dstNode->addArc(a_srcBlock, dir_i_10, dir_i_01);
    }
}

unsigned int MBGraph::numBoundaries(unsigned int a_index, Point& a_dir) const
{
    PROTO_ASSERT(a_index < size(),
            "MBGraph::numBoundaries | Error: Block index %u is out of bounds (size = %u).",
            a_index, size());
    auto arcs = m_blocks[a_index]->arcs;
    if (arcs.find(a_dir) == arcs.end()) { return 0; }
    return m_blocks[a_index]->arcs[a_dir].size();
}

std::set<Point> MBGraph::fullConnectivity(unsigned int a_srcBlock, unsigned int a_dstBlock) const
{
    auto& srcNode = *m_blocks[a_srcBlock];
    std::set<Point> ret;
    Box K = Box::Kernel(1);
    for (auto dir : K)
    {
        if (dir == Point::Zeros()) { continue; }
        for (auto arc : srcNode.arcs[dir])
        {
            if (arc.dstBlock == a_dstBlock)
            {
                ret.insert(dir);
            }
        }
    }
    return ret; 
}

Point MBGraph::connectivity(unsigned int a_srcBlock, unsigned int a_dstBlock) const
{
    if (a_srcBlock == a_dstBlock) { return Point::Zeros(); }
    auto conns = fullConnectivity(a_srcBlock, a_dstBlock);
    unsigned int minCodim = DIM+1;
    Point minDir = Point::Zeros();
    for (auto dir : conns)
    {
        int thisCodim = dir.codim();
        if (thisCodim < minCodim)
        {
            minDir = dir;
            minCodim = thisCodim;
        }
    }
    return minDir;
}

Point MBGraph::reverseArc(
        unsigned int a_srcBlock,
        unsigned int a_dstBlock,
        Point a_dir) const
{
    PROTO_ASSERT(a_srcBlock < size(),
        "MBGraph::reverseArc | Error: Block index %u is out of bounds.", a_srcBlock);
    PROTO_ASSERT(a_dstBlock < size(),
        "MBGraph::reverseArc | Error: Block index %u is out of bounds.", a_dstBlock);
    auto& bounds = m_blocks[a_srcBlock]->arcs;
    PROTO_ASSERT(bounds.find(a_dir) != bounds.end(),
        "MBGraph::reverseArc | Error: Block %u has no boundary in the specified direction",
        a_srcBlock);
    auto& arcs = bounds[a_dir];
    for (auto arc : arcs)
    {
        if (arc.dstBlock == a_dstBlock) { return arc.dstToSrc; }
    }
    return Point::Zeros();
}

Point MBGraph::reverseDir(
        unsigned int a_srcBlock,
        unsigned int a_dstBlock,
        Point a_dir) const
{
    PROTO_ASSERT(a_srcBlock < size(),
        "MBGraph::reverseDir | Error: Block index %u is out of bounds.", a_srcBlock);
    PROTO_ASSERT(a_dstBlock < size(),
        "MBGraph::reverseDir | Error: Block index %u is out of bounds.", a_dstBlock);
    auto bounds = m_blocks[a_srcBlock]->arcs;
    Point norm = connectivity(a_srcBlock, a_dstBlock);
    auto arcs = bounds[norm];
    Point adjNorm = Point::Ones();
    for (auto arc : arcs)
    {
        if (arc.dstBlock == a_dstBlock)
        {
            adjNorm = arc.dstToSrc;
        }
    }
    
    if (a_dir == norm) { return adjNorm; }

    auto R = rotation(a_srcBlock, norm, a_dstBlock);
    Point orth;
    for (int nn = 0; nn < DIM; nn++)
    {
        orth[nn] = (norm[nn] == 0) ? 1 : 0;
    }
    return -R(a_dir*orth) + adjNorm;
}

inline CoordPermutation& MBGraph::rotation(unsigned int a_block, Point a_dir) const
{
    PROTO_ASSERT(a_block < size(),
            "MBGraph::rotation | Error: Block index %u is out of bounds.", a_block);
    PROTO_ASSERT(a_dir.codim() == 1,
            "MBGraph::rotation | Error: Directional input must be a codimension 1 boundary.");
    PROTO_ASSERT(numBoundaries(a_block, a_dir) == 1,
            "MBGraph::rotation | Error: Block %u has no boundary in the specified direction.",
            a_block);
    return m_blocks[a_block]->arcs[a_dir][0].R;
}

/// START HERE
inline CoordPermutation& MBGraph::rotation(
        unsigned int a_srcBlock,
        Point a_dir,
        unsigned int a_dstBlock) const
{
    PROTO_ASSERT(a_srcBlock < size(),
            "MBGraph::rotation | Error: Block index %u is out of bounds.", a_srcBlock);
    PROTO_ASSERT(a_dstBlock < size(),
            "MBGraph::rotation | Error: Block index %u is out of bounds.", a_dstBlock);
    PROTO_ASSERT(a_dir.codim() > 0,
            "MBGraph::rotation | Error: Directional input has invalid codimension.");
    PROTO_ASSERT(numBoundaries(a_srcBlock, a_dir) > 0,
            "MBGraph::rotation | Error: Block %u has no boundary in the specified direction.",
            a_srcBlock);
    if (a_dir.codim() == 1)
    {
        return rotation(a_srcBlock, a_dir);
    } else {
        return m_blocks[a_srcBlock]->rotation(a_dir, a_dstBlock);
    }
}

unsigned int MBGraph::adjacent(
        unsigned int    a_srcBlock,
        unsigned int    a_dir,
        Side::LoHiSide  a_side) const
{
    return adjacent(a_srcBlock, Point::Basis(a_dir, a_side));
}

unsigned int MBGraph::adjacent(
        unsigned int    a_srcBlock,
        Point           a_dir) const
{
    PROTO_ASSERT(a_dir.codim() == 1,
            "MBGraph::adjacent | Error: Directional input must have codimension 1.");
    if (numBoundaries(a_srcBlock, a_dir) != 1) { return size(); }
    return m_blocks[a_srcBlock]->arcs[a_dir][0].dstBlock;
}

std::vector<MBGraphArc> MBGraph::boundaries(unsigned int a_block, Point a_dir) const
{
    auto arcs = m_blocks[a_block]->arcs[a_dir];
    std::vector<MBGraphArc> bounds;
    for (auto arc : arcs)
    {
        unsigned int dstBlock = arc.dstBlock;
        if (connectivity(a_block, dstBlock).codim() == a_dir.codim())
        {
            bounds.push_back(arc);
        }
    }
    return bounds;
}

bool MBGraph::operator==(const MBGraph& a_rhs) const
{
    // TODO: Most applications won't need more than one graph
    //  We can change this in the future if need be. -CLG
    return (this == &a_rhs);
}

void MBGraph::print() const
{
    using Proto::pout;
    for (auto block : m_blocks)
    {
        pout() << "Block " << block->index << std::endl;
        Box K = Box::Kernel(1);
        for (auto dir : K)
        {
            pout() << "\tDir: " << dir << " (codim = " << dir.codim() << ")" << std::endl;
            for (auto arc : block->arcs[dir])
            {
                pout() << "\t\tBound: Block " << arc.dstBlock << " | Local Dir: " << arc.srcToDst << std::endl;
            }
        }
    }
}

void MBGraph::close()
{
    if (m_closed) { return; }
    fixRotations();
    for (int bi = 0; bi < size(); bi++)
    {
        for (auto dir_ij : Box::Kernel(1))
        {
            if (dir_ij == Point::Zeros()) { continue; }
            for (auto bound_ij : boundaries(bi, dir_ij))
            {
                PROTO_ASSERT(bound_ij.srcToDst == dir_ij,
                        "MBGraph::close() | Error: Graph corruption");
                PROTO_ASSERT(bound_ij.srcBlock == bi,
                        "MBGraph::close() | Error: Graph corruption");
                int bj = bound_ij.dstBlock;
                auto dir_ji = bound_ij.dstToSrc;
                if (dir_ij.codim() == 1)
                {
                    PROTO_ASSERT(adjacent(bi,dir_ij) == bj,
                        "MBGraph::close() | Error: Graph corruption");
                    PROTO_ASSERT(adjacent(bj,dir_ji) == bi,
                        "MBGraph::close() | Error: Graph corruption");
                }
                if (dir_ij.codim() <= 2)
                {
                    auto R = rotation(bi, dir_ij, bj);
                    PROTO_ASSERT(R(dir_ij) == -dir_ji,
                        "MBGraph::close() | Error: Graph corruption");
                }
            }
        }
    }
}




