
template<typename T, unsigned int C, MemType MEM, Centering CTR>
MBLevelBoxData<T, C, MEM, CTR>::MBLevelBoxData(
        MBDisjointBoxLayout& a_layout,
        Point a_ghost)
{
    std::array<Point, DIM+1> ghost;
    ghost.fill(a_ghost);
    define(a_layout, ghost);
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
MBLevelBoxData<T, C, MEM, CTR>::MBLevelBoxData(
        MBDisjointBoxLayout&     a_layout,
        std::array<Point, DIM+1> a_ghost)
{
    define(a_layout, a_ghost);
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void MBLevelBoxData<T, C, MEM, CTR>::define(
        MBDisjointBoxLayout&      a_layout,
        std::array<Point, DIM+1> a_ghost)
{
    // a_ghost convention:
    // [interior, codim1, codim2, ...]
    m_ghost = a_ghost;
    m_layout = &a_layout;
    m_data.clear();
    m_data.resize(a_layout.numBlocks());
    for (int bi = 0; bi < a_layout.numBlocks(); bi++)
    {
        m_data[bi] = std::make_shared<LevelBoxData<T, C, MEM, CTR>>(
                a_layout.blockLayout(bi), a_ghost[0]);
    }
    for (auto iter : a_layout)
    {
        auto block = a_layout.block(iter);
        auto patch = a_layout.point(iter);
        auto localLayout = a_layout.blockLayout(block);
        auto patchDomain = localLayout.patchDomain().box();
        auto localIndex = localLayout.find(patch);
        auto patchBox = localLayout[localIndex];
        

        Box K = Box::Kernel(1);
        for (int cc = 1; cc <= DIM; cc++)
        {
            int ghostSize = m_ghost[cc].max();
            for (auto di : K)
            {
                if (codim(di) != cc) { continue; }
                if (patchDomain.contains(patch+di)) { continue; }
                Point shrinkDirs;
                for (int d = 0; d < DIM; d++)
                {
                    shrinkDirs[d] = (di[d] != 0) ? -1 : 0;
                }
                Box k = K.grow(shrinkDirs).shift(di);
                PROTO_ASSERT(k.contains(di),
                        "MBLevelBoxData::define | Error: Data corruption");
                auto bounds = a_layout.domain().boundaries(block, di);
                Box boundaryPatches = patchDomain.adjacent(di,1);
                for (auto bound : bounds)
                {
                    PROTO_ASSERT(bound.srcToDst == di,
                        "MBLevelBoxData::define | Error: Data corruption");
                    auto adjBlock = bound.dstBlock;
                    Point adjDir = bound.dstToSrc;
                    auto adjLocalLayout = a_layout.blockLayout(adjBlock);
                    Box adjPatchDomain = adjLocalLayout.patchDomain().box();
                    Box adjBoundaryPatches = adjPatchDomain.edge(adjDir);
                    auto R = a_layout.domain().graph().rotation(block, di, adjBlock);
                    for (auto ki : k)
                    {
                        boundsKey_t key(patch, block, ki);
                        Point neighbor = patch + ki;
                        if (boundaryPatches.contains(neighbor))
                        {
                            Point L = boundaryPatches.low();
                            auto adjPatch = R.rotateCell(neighbor - L, boundaryPatches.shift(-L))
                                + adjBoundaryPatches.low();
                            auto adjIndex = adjLocalLayout.find(adjPatch);
                            
                            // for codim > 1, the rotation isn't defined in the directions
                            // that are linearly independent from the boundary surface
                            // to get this right, we rotate the in-plane coordinates and 
                            // use the reverse arc for the out of plane coordinates
                            Point norm = a_layout.domain().graph().connectivity(block, adjBlock);
                            Point adjNorm = a_layout.domain().graph().reverseArc(block, adjBlock, norm);
                            Point orth;
                            for (int nn = 0; nn < DIM; nn++)
                            {
                                orth[nn] = (norm[nn] == 0) ? 1 : 0;
                            }
                            Point kiAdj = -R(ki*orth) + adjNorm;

                            //TODO: make this work with anisotropic ghost regions
                            Box patchBoundary = patchBox.adjacent(ki, ghostSize);
                            //Box patchBoundary = patchBox.adjacent(ki * m_ghost[cc]);
                            Box adjPatchBoundary = adjLocalLayout[adjIndex].edge(kiAdj,ghostSize);
                            //Box adjPatchBoundary = adjLocalLayout[adjIndex].edge(kiAdj * m_ghost[cc]);
                             
                            auto adjIter = a_layout.find(adjPatch, adjBlock);
                            PROTO_ASSERT(adjIter != *a_layout.end(),
                                    "MBLevelBoxData::define | Error: Attempting to create \
                                    BoundaryData with nonexistent patch");
                            MBBoundaryData<T, C, MEM> boundData(
                                    iter, adjIter, patchBoundary, adjPatchBoundary);
                            m_bounds[key].push_back(boundData);
                        }
                    }
                }
            }
        }
    }
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
template<typename Func, typename... Srcs>
void MBLevelBoxData<T, C, MEM, CTR>::initialize(Func& a_func, Srcs... a_srcs)
{
    for (int bi = 0; bi < numBlocks(); bi++)
    {
        auto& data = blockData(bi);
        data.initialize(a_func, bi, a_srcs...);
    }
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
template<typename Func, typename... Srcs>
void MBLevelBoxData<T, C, MEM, CTR>::initConvolve(Func& a_func, Srcs... a_srcs)
{
    for (int bi = 0; bi < numBlocks(); bi++)
    {
        auto& data = blockData(bi);
        data.initConvolve(a_func, bi, a_srcs...);
    }
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void
MBLevelBoxData<T, C, MEM, CTR>::fillBoundaries()
{

}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void
MBLevelBoxData<T, C, MEM, CTR>::interpolate()
{
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void
MBLevelBoxData<T, C, MEM, CTR>::exchange()
{
    fillBoundaries();
    interpolate();
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
BoxData<T, C, MEM>&
MBLevelBoxData<T, C, MEM, CTR>::operator[](MBIndex a_index)
{
    unsigned int block = m_layout->block(a_index);
    auto localIndex = m_layout->blockIndex(a_index);
    return (*m_data[block])[localIndex];
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
const BoxData<T, C, MEM>&
MBLevelBoxData<T, C, MEM, CTR>::operator[](MBIndex a_index) const
{
    unsigned int block = m_layout->block(a_index);
    auto localIndex = m_layout->blockIndex(a_index);
    return (*m_data[block])[localIndex];
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
std::array<Point, DIM+1> 
MBLevelBoxData<T, C, MEM, CTR>::ghost() const
{
    return m_ghost;
}
template<typename T, unsigned int C, MemType MEM, Centering CTR>
LevelBoxData<T, C, MEM, CTR>&
MBLevelBoxData<T, C, MEM, CTR>::blockData(unsigned int a_block)
{
    return *m_data[a_block];
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
const LevelBoxData<T, C, MEM, CTR>&
MBLevelBoxData<T, C, MEM, CTR>::blockData(unsigned int a_block) const
{
    return *m_data[a_block];
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
std::vector<MBBoundaryData<T, C, MEM>>&
MBLevelBoxData<T, C, MEM, CTR>::bounds(MBIndex a_index, Point a_dir)
{
    auto patch = m_layout->point(a_index);
    auto block = m_layout->block(a_index);
    return m_bounds[boundsKey_t(patch, block, a_dir)];
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
const std::vector<MBBoundaryData<T, C, MEM>>&
MBLevelBoxData<T, C, MEM, CTR>::bounds(MBIndex a_index, Point a_dir) const
{
    auto patch = m_layout->point(a_index);
    auto block = m_layout->block(a_index);
    return m_bounds[boundsKey_t(patch, block, a_dir)];
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
MBIterator
MBLevelBoxData<T, C, MEM, CTR>::begin()
{
    return m_layout->begin();
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
MBIterator
MBLevelBoxData<T, C, MEM, CTR>::end()
{
    return m_layout->end();
}

//=======================================================================================
/// LEVEL EXCHANGE COPIER OP

template<typename T, unsigned int C, MemType MEM, Centering CTR>
MBLevelExchangeCopierOp<T, C, MEM, CTR>::MBLevelExchangeCopierOp(
        MBLevelBoxData<T, C, MEM, CTR>& a_data)
{
    m_data = &a_data;
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
uint64_t
MBLevelExchangeCopierOp<T, C, MEM, CTR>::linearSize(
        const Box& a_box,
        const MBIndex& a_index) const
{
    //TODO: Does this do the right thing?
    return C*a_box.size()*sizeof(T);
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void
MBLevelExchangeCopierOp<T, C, MEM, CTR>::linearOut(
        void*           a_buf, 
        const Box&      a_bx,
        const MBIndex   a_index) const
{
    MemType mem = pointerMemType(a_buf);
    PROTO_ASSERT(mem == MEM,
        "MBLevelExchangeCopierOp::linearOut | Error: Memtype mismatch. Detected %i but should be %i",
        mem, MEM);
    auto& patch = m_data->operator[](a_index);
    CInterval cint(0,C-1);
    patch.linearOut(a_buf, a_bx, cint);
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void
MBLevelExchangeCopierOp<T, C, MEM, CTR>::linearIn(
        void*           a_buf, 
        const Box&      a_bx,
        const MBIndex   a_index) const
{
    MemType mem = pointerMemType(a_buf);
    PROTO_ASSERT(mem == MEM,
        "MBLevelCopierOp::linearIn | Error: Memtype mismatch. Detected %i but should be %i",
        mem, MEM);
    Box dirs = Box::Kernel(1);
    CInterval cint(0,C-1);
    for (auto dir : dirs)
    {
        auto bounds = m_data->bounds(a_index, dir);
        for (auto bound : bounds)
        {
            auto& boundData = *bound.localData;
            if (boundData.box() == a_bx)
            {
                boundData.linearIn(a_buf, a_bx, cint);
                return;
            }
        }
    }
    MayDay<void>::Warning("MBLevelExchangeCopierOp::linearIn | Warning: \
            could not find a boundary to copy data into");
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void
MBLevelExchangeCopierOp<T, C, MEM, CTR>::localCopy(
        const Box&      a_domain,
        const MBIndex   a_domainIndex,
        const Box&      a_range, 
        const MBIndex   a_rangeIndex) const
{
    bool linearizationTest = false;
    if (!linearizationTest)
    {
        /*
        PR_TIMERS("MBLevelCopierOp::localCopy");
        auto& src = m_data->operator[](a_domainIndex);
        auto srcBlock = m_data->layout().block(a_domainIndex);
        Box dirs = Box::Kernel(1);
        for (auto dir : dirs)
        {
            auto& bound = m_data->bound(a_rangeIndex, dir, srcBlock);
            Box b = bound.adjData->box() & a_domain
            if (bound.adjData->box() & a_domain))
            {

            }
        }
        auto& dst = m_data->operator[](a_rangeIndex);
        Point shift = a_range.low() - a_domain.low();
        src.copyTo(dst, a_domain, shift);
        (*/
    }
    else
    {
        PR_TIMERS("MBLevelCopierOp::localCopy (linearization test)");
        int sizeSrc = linearSize(a_domain, a_domainIndex);
        int sizeDst = linearSize(a_range, a_rangeIndex);

        if (sizeSrc != sizeDst)
        {
            MayDay<void>::Abort(
                    "LinearizationTest failure: dest and source have different sizes");
        }
        
        void* srcBuffer = proto_malloc<MEM>(sizeSrc);
        void* dstBuffer = proto_malloc<MEM>(sizeDst);
        {
            PR_TIMERS("linearIn / linearOut inside LDOP::op");
            linearOut(srcBuffer, a_domain, a_domainIndex);
            proto_memcpy<MEM, MEM>(srcBuffer, dstBuffer, sizeSrc);
            linearIn (dstBuffer, a_range,  a_rangeIndex);
        }
    }
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void
MBLevelExchangeCopier<T, C, MEM, CTR>::buildMotionPlans(
        MBLevelExchangeCopierOp<T, C, MEM, CTR>& a_op)
{
    auto layout = a_op.m_src->layout();
    auto ghost  = a_op.m_src->ghost();

    Box dirs = Box::Kernel(1);
    int localRank = procID();
    for (auto iter : layout)
    {
        for (auto dir : dirs)
        {
            auto bounds = a_op.m_src.bounds(iter, dir);
            for (auto bound : bounds)
            {
                PROTO_ASSERT(layout.procID(bound.localIndex) == localRank,
                        "MBLevelExchangeCopier::buildMotionPlans | Error: Data corruption.");
                int adjRank = layout.procID(bound.adjIndex);
                Box adjBox = bound.adjData->box();
                Box localBox = bound.localData->box();
                MotionItem<MBBoxPartition> toMotionItem(
                        bound.adjIndex, bound.localIndex, adjBox, localBox);
                if (localRank == adjRank)
                {
                    this->m_localMotionPlan.push_back(toMotionItem);
                } else {
                    toMotionItem.procID = adjRank;
                    this->m_toMotionPlan.push_back(toMotionItem);
                    MotionItem<MBBoxPartition> fromMotionItem(
                            bound.localIndex, bound.adjIndex, localBox, adjBox);
                    fromMotionItem.procID = adjRank;
                    this->m_fromMotionPlan.push_back(fromMotionItem);
                }
            }
        }
    }
    this->sort();
    /*
    pout() << std::endl << "=======================================================================" << std::endl;
    pout() << "Printing ExchangeCopier" << std::endl;
    pout() << "Layout" << std::endl;
    layout.print();
    pout() << "Copier" << std::endl;
    this->print();
    */
}

template<typename T, unsigned int C, MemType SRC_MEM, MemType DST_MEM, Centering CTR>
MBLevelCopierOp<T, C, SRC_MEM, DST_MEM, CTR>::MBLevelCopierOp(
        MBLevelBoxData<T, C, SRC_MEM, CTR>& a_src,
        MBLevelBoxData<T, C, DST_MEM, CTR>& a_dst)
{
    PROTO_ASSERT(SRC_MEM == DST_MEM,
        "MBLevelCopierOp::constructor | Error: \
        MBLevelCopierOp with different MemType src and dst is not yet implemented.");
    PROTO_ASSERT(a_src.layout().domain()==a_dst.layout().domain(),
            "MBLevelCopier::Constructor | Error: \
            src and dst must have the same problem domain");
    m_src = &a_src;
    m_dst = &a_dst;
}

template<typename T, unsigned int C, MemType SRC_MEM, MemType DST_MEM, Centering CTR>
uint64_t
MBLevelCopierOp<T, C, SRC_MEM, DST_MEM, CTR>::linearSize(
        const Box& a_box,
        const MBIndex a_index) const
{
    //TODO: Does this do the right thing?
    return C*a_box.size()*sizeof(T);
}

template<typename T, unsigned int C, MemType SRC_MEM, MemType DST_MEM, Centering CTR>
void
MBLevelCopierOp<T, C, SRC_MEM, DST_MEM, CTR>::linearOut(
        void*           a_buf, 
        const Box&      a_bx,
        const MBIndex  a_index) const
{
    MemType mem = pointerMemType(a_buf);
    PROTO_ASSERT(mem == SRC_MEM,
        "MBLevelCopierOp::linearOut | Error: Memtype mismatch. Detected %i but should be %i",
        mem, SRC_MEM);
    auto& patch = m_src->operator[](a_index);
    CInterval cint(0,C-1);
    patch.linearOut(a_buf, a_bx, cint);
}

template<typename T, unsigned int C, MemType SRC_MEM, MemType DST_MEM, Centering CTR>
void
MBLevelCopierOp<T, C, SRC_MEM, DST_MEM, CTR>::linearIn(
        void*           a_buf, 
        const Box&      a_bx,
        const MBIndex  a_index) const
{
    MemType mem = pointerMemType(a_buf);
    PROTO_ASSERT(mem == DST_MEM,
        "MBLevelCopierOp::linearIn | Error: Memtype mismatch. Detected %i but should be %i",
        mem, DST_MEM);
    auto& patch = m_dst->operator[](a_index);
    CInterval cint(0,C-1);
    patch.linearIn(a_buf, a_bx, cint);
}

template<typename T, unsigned int C, MemType SRC_MEM, MemType DST_MEM, Centering CTR>
void
MBLevelCopierOp<T, C, SRC_MEM, DST_MEM, CTR>::localCopy(
        const Box& a_domain,
        const MBIndex a_domainIndex,
        const Box& a_range, 
        const MBIndex a_rangeIndex) const
{
    bool linearizationTest = false;
    if (!linearizationTest)
    {
        PR_TIMERS("MBLevelCopierOp::localCopy");
        auto& src = m_src->operator[](a_domainIndex);
        auto& dst = m_dst->operator[](a_rangeIndex);
        Point shift = a_range.low() - a_domain.low();
        src.copyTo(dst, a_domain, shift);
    }
    else
    {
        PR_TIMERS("MBLevelCopierOp::localCopy (linearization test)");
        int sizeSrc = linearSize(a_domain, a_domainIndex);
        int sizeDst = linearSize(a_range, a_rangeIndex);

        if (sizeSrc != sizeDst)
        {
            MayDay<void>::Abort(
                    "LinearizationTest failure: dest and source have different sizes");
        }
        
        void* srcBuffer = proto_malloc<SRC_MEM>(sizeSrc);
        void* dstBuffer = proto_malloc<DST_MEM>(sizeDst);
        {
            PR_TIMERS("linearIn / linearOut inside LDOP::op");
            linearOut(srcBuffer, a_domain, a_domainIndex);
            proto_memcpy<SRC_MEM, DST_MEM>(srcBuffer, dstBuffer, sizeSrc);
            if (SRC_MEM != DST_MEM)
            {
#ifdef PROTO_CUDA
                protoThreadSynchronize();
#endif
            }
            linearIn (dstBuffer, a_range,  a_rangeIndex);
        }
    }
}


/// LEVEL COPIER (WORK IN PROGRESS)
template<typename T, unsigned int C, MemType SRC_MEM, MemType DST_MEM, Centering CTR>
void
MBLevelCopier<T, C, SRC_MEM, DST_MEM, CTR>::buildMotionPlans(MBLevelCopierOp<T, C, SRC_MEM, DST_MEM, CTR>& a_op)
{
    PROTO_ASSERT(SRC_MEM == DST_MEM, 
        "MBLevelCopier::buildMotionPlans | Error: \
        Copy operations between HOST and DEVICE are not yet implemented.");
    // Construct Motion items for copyTo and insert into MotionPlans.
    // For each motion item, the source Box is specified relative to the location
    // of the source patch, and the dest Box relative to the location of the dest
    // patch. This is easy for the fromMotionPlan / localMotionPlan. In constructing
    // the toMotionPlan, though, must explicitly shift back and forth between local
    // coordinate systems defined by the source and destination.
    int destProcID = procID();

    auto srclayout = a_op.m_src->layout();
    auto destlayout = a_op.m_dst->layout();
    auto domain = srclayout.domain();
    auto destGhost = a_op.m_dst->ghost();
    // For each Box on this processor, compute all the pairs of motionItems that move data
    // to this Box.
    for (auto diter : destlayout) 
    {
        // destGhost[0] is the ghost size in the bulk of each block
        Box bx = destlayout[diter].grow(destGhost[0]);

        // Compute all the patches (as Points in a bitmap) that might move data to bx.

        Box srcPatches = bx.coarsen(srclayout.boxSize());
        // Iterate over the possible source patches to contruct MotionItems.

        for (auto pt : srcPatches)
        {
            // Find the possible source patch, or its periodic image.
            
            auto di = srclayout.find(pt);
            if (di != *srclayout.end())
            {
                Box destBx = Box(pt, pt).refine(srclayout.boxSize()) & bx;
                // Compute srcBx to be the shifted image of destBx. If destBx is in
                // the ProblemDomain box, srcBx == destBx.
                //Box srcBx = Box(domain.image(destBx.low()),domain.image(destBx.high()));

                if (!destBx.empty())
                {
                    // Create MotionItem.
                    // Compute srcBx to be the shifted image of destBx. If destBx 
                    // is in the ProblemDomain box, srcBx == destBx.
                    Box srcBx = Box(domain.image(destBx.low()),
                                    domain.image(destBx.high()));
                    MotionItem<BoxPartition> toMotionItem(di, diter, srcBx, destBx);
                    int srcProcID = srclayout.procID(di);
                    if (srcProcID == procID())
                    {
                        PROTO_ASSERT(srcProcID == destlayout.procID(diter),
                            "LevelCopier::buildMotionPlan | Error: Data corruption.");
                        this->m_localMotionPlan.push_back(toMotionItem);
                    }
                    else
                    {
                        PROTO_ASSERT(srcProcID != destlayout.procID(diter),
                            "LevelCopier::buildMotionPlan | Error: Data corruption.");
                        toMotionItem.procID = srcProcID;
                        this->m_toMotionPlan.push_back(toMotionItem);
                    }
                }
            }
        }
    }

    // Iterate over Boxes in the source to construct fromMotionPlan  
    for (auto siter : srclayout)
    {
        // We grow the source Box by dstGhost to capture all the patches that the
        // source will copy to. dProbC contains all the points in the bitmap for
        // the destination.
        Box srcBox = srclayout[siter];
        Box dProbC = srcBox.grow(destGhost).coarsen(destlayout.boxSize());
        for (auto pt : dProbC)
        {
            auto di = destlayout.find(pt);
            if (di != *destlayout.end())
            {
                Box bxDestShifted = Box(pt, pt).refine(destlayout.boxSize()).grow(destGhost);
                Box bxDest = destlayout[di].grow(destGhost);
                Point shiftSrcDest = bxDest.low() - bxDestShifted.low();
                bxDest = bxDest & (srcBox.shift(shiftSrcDest));
                Box bxSrc = bxDestShifted & srcBox;
                if (!bxSrc.empty())
                {
                    int destProcID = destlayout.procID(di);           
                    if (destProcID != procID())
                    {
                        MotionItem<BoxPartition> fromMotionItem(siter, di, bxSrc, bxDest);
                        fromMotionItem.procID =destProcID;
                        this->m_fromMotionPlan.push_back(fromMotionItem);
                    }
                }
            }
        }
    } // end src loop
    this->sort();
    /*
    pout() << std::endl << "=======================================================================" << std::endl;
    pout() << "Printing LevelCopier" << std::endl;
    pout() << "Source Layout" << std::endl;
    srclayout.print();
    pout() << "Dest Layout" << std::endl;
    destlayout.print();
    this->print();
    */
}

