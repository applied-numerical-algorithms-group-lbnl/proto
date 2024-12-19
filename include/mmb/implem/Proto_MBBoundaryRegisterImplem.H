
template<typename T, unsigned int C, MemType MEM, Centering CTR>
MBBoundaryRegister<T,C,MEM,CTR>::MBBoundaryRegister(
        MBDisjointBoxLayout         a_layout,   // layout
        int                         a_depth,    // distance from boundary
        Point                       a_ghost,    
        bool                        a_bothSides)
{
    define(a_layout, a_depth, a_ghost, a_bothSides);
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void
MBBoundaryRegister<T,C,MEM,CTR>::define(
        MBDisjointBoxLayout         a_layout,   // layout
        int                         a_depth,    // distance from boundary
        Point                       a_ghost,    // additional space around each register
        bool                        a_bothSides)// allocate space on both "sides" of the boundary
{
    m_layout = a_layout;
    m_bounds.clear();

    for (auto iter : a_layout)
    {
        auto block = a_layout.block(iter);
        auto patch = a_layout.point(iter);
        auto localLayout = a_layout.getBlock(block);
        auto patchDomain = localLayout.patchDomain().box();
        auto localIndex = localLayout.find(patch);
        auto patchBox = localLayout[localIndex];

        Box K = Box::Kernel(1);
        for (int cc = 1; cc <= DIM; cc++)
        {
            for (auto di : K)
            {
                // skip this codimension
                if (di.codim() != cc)
                {
                    continue;
                }
                // adjacent patch is in the same block
                if (patchDomain.containsPoint(patch + di))
                {
                    continue;
                }
                // all patches that could be on the boundary of this one (including on diag)
                Point shrinkDirs;
                for (int d = 0; d < DIM; d++)
                {
                    shrinkDirs[d] = (di[d] != 0) ? -1 : 0;
                }
                Box k = K.grow(shrinkDirs).shift(di);
                PROTO_ASSERT(k.containsPoint(di),
                             "MBLevelBoxData::define | Error: Data corruption");
                // all patches that bound this block in the given direction
                Box boundaryPatches = patchDomain.adjacent(di, 1);
                
                // loop over all blocks that bound this one in direction di
                auto bounds = a_layout.domain().boundaries(block, di);
                for (auto bound : bounds)
                {
                    PROTO_ASSERT(bound.srcToDst == di,
                                 "MBLevelBoxData::define | Error: Data corruption");
                    auto adjBlock = bound.dstBlock;
                    Point adjDir = bound.dstToSrc;
                    auto& adjLocalLayout = a_layout.getBlock(adjBlock);
                    Box adjPatchDomain = adjLocalLayout.patchDomain().box();
                    Box adjBoundaryPatches = adjPatchDomain.edge(adjDir);
                    auto R = a_layout.domain().graph().rotation(block, di, adjBlock);
                    for (auto ki : k)
                    {
                        boundsKey_t key(patch, block, ki);
                        Point neighbor = patch + ki;
                        // found a valid potential neighbor (if it exists)
                        if (boundaryPatches.containsPoint(neighbor))
                        {
                            auto adjPatch = R.rotateCell(
                                    neighbor, boundaryPatches, adjBoundaryPatches);
                            auto adjIndex = adjLocalLayout.find(adjPatch);

                            // check if the adjacent layout contains the desired patch
                            if (adjIndex == adjLocalLayout.end()) { continue; }
                            Box adjPatchBox = adjLocalLayout[adjIndex];
                            // for codim > 1, the rotation isn't defined in the directions
                            // that are linearly independent from the boundary surface
                            // to get this right, we rotate the in-plane coordinates and
                            // use the reverse arc for the out of plane coordinates
                            
                            Point kiAdj = a_layout.domain().graph().reverseDir(
                                    block, adjBlock, ki);
                            Box patchBoundary = patchBox.adjacent(ki, a_depth);
                            Box adjPatchBoundary;
                            if (a_bothSides)
                            {
                                int sign = a_depth / abs(a_depth);
                                patchBoundary = patchBoundary.extrude(ki*(-sign),abs(a_depth)); 
                                adjPatchBoundary = patchBoundary; 
                            } else {
                                adjPatchBoundary = patchBox.edge(ki, a_depth);
                            }

                            PROTO_ASSERT(patchBoundary.size() == adjPatchBoundary.size(), 
                                    "MBBoundaryRegister::define | Error: data corruption.");

                            auto adjIter = a_layout.find(adjPatch, adjBlock);
                            
                            // redundant check if the adjacent layout contains the desired patch
                            if (adjIter == *a_layout.end()) { continue; }
                            MBBoundaryData<T, C, MEM> boundData(a_layout,
                                iter, adjIter, patchBoundary, adjPatchBoundary,
                                R.inverse(), a_ghost);
                            m_bounds[key].push_back(boundData);
                        }
                    }
                }
            }
        }
    }
    m_copier.define(MBBoundaryCopierOp<T,C,MEM,CTR>(*this));
}

// Get All Boundary Data
template<typename T, unsigned int C, MemType MEM, Centering CTR>
std::vector<MBBoundaryData<T,C,MEM>> MBBoundaryRegister<T,C,MEM,CTR>::bounds(
        MBIndex a_index,
        Point   a_dir)
{
    std::vector<MBBoundaryData<T,C,MEM>> empty;
    auto block = m_layout.block(a_index);
    auto patch = m_layout.point(a_index);
    boundsKey_t key(patch, block, a_dir);
    if (m_bounds.find(key) == m_bounds.end())
    {
        return empty;
    }
    return m_bounds[key];
}

// Get All Boundary Data
template<typename T, unsigned int C, MemType MEM, Centering CTR>
std::vector<MBBoundaryData<T,C,MEM>>
MBBoundaryRegister<T,C,MEM,CTR>::bounds(
        MBIndex a_index)
{
    std::vector<MBBoundaryData<T,C,MEM>> ret;
    for (auto dir : Box::Kernel(1))
    {
        auto bi = bounds(a_index, dir);
        ret.insert(ret.end(), bi.begin(), bi.end());
    }
    return ret;
}

// Get Boundary Data
template<typename T, unsigned int C, MemType MEM, Centering CTR>
MBBoundaryData<T,C,MEM>&
MBBoundaryRegister<T,C,MEM,CTR>::bound(
        MBIndex         a_index,
        Point           a_dir,
        unsigned int    a_adjBlock)
{
    for (auto bound : bounds(a_index, a_dir))
    {
        if (bound.dstBlock == a_adjBlock)
        {
            return bound;
        }
    }
    MayDay<void>::Error("MBBoundaryRegister::bound | boundary not found.");
}

// Get Boundary Data
template<typename T, unsigned int C, MemType MEM, Centering CTR>
MBBoundaryData<T,C,MEM>& MBBoundaryRegister<T,C,MEM,CTR>::bound(
        MBIndex         a_index,
        unsigned int    a_dir,
        Side::LoHiSide  a_side)
{
    auto boundaries = bounds(a_index, Point::Basis(a_dir, a_side));
    PROTO_ASSERT(boundaries.size() == 1,
            "MBBoundaryRegister::bound | Error: boundary not found.");
    return boundaries[0];
}

// Get Local Data
template<typename T, unsigned int C, MemType MEM, Centering CTR>
BoxData<T,C,MEM>& MBBoundaryRegister<T,C,MEM,CTR>::local(
        MBIndex         a_index,
        unsigned int    a_dir,
        Side::LoHiSide  a_side)
{
    return *(bound(a_index, a_dir, a_side).localData);
}

// Get Local Data
template<typename T, unsigned int C, MemType MEM, Centering CTR>
BoxData<T,C,MEM>& MBBoundaryRegister<T,C,MEM,CTR>::local(
        MBIndex         a_index,
        Point           a_dir,
        unsigned int    a_adjBlock)
{
    return *(bound(a_index, a_dir, a_adjBlock).localData);
}

// Get Local Data
template<typename T, unsigned int C, MemType MEM, Centering CTR>
BoxData<T,C,MEM>& MBBoundaryRegister<T,C,MEM,CTR>::local(
        MBIndex         a_localIndex,
        MBIndex         a_adjIndex)
{
    Box dirs = Box::Kernel(1);
    for (auto dir : dirs)
    {
        for (auto bound : bounds(a_localIndex, dir))
        {
            if (bound.adjIndex == a_adjIndex)
            {
                return *(bound.localData);
            }
        }
    }
    MayDay<void>::Error("MBBoundaryRegister::local | Error: no boundary between patches.");
}

// Get Adjacent Data
template<typename T, unsigned int C, MemType MEM, Centering CTR>
BoxData<T,C,MEM>& MBBoundaryRegister<T,C,MEM,CTR>::adjacent(
        MBIndex         a_index,
        unsigned int    a_dir,
        Side::LoHiSide  a_side)
{
    return *(bound(a_index, a_dir, a_side).adjData);
}

// Get Adjacent Data
template<typename T, unsigned int C, MemType MEM, Centering CTR>
BoxData<T,C,MEM>& MBBoundaryRegister<T,C,MEM,CTR>::adjacent(
        MBIndex         a_index,
        Point           a_dir,
        unsigned int    a_adjBlock)
{
    return *(bound(a_index, a_dir, a_adjBlock).adjData);
}

// Get Adjacent
template<typename T, unsigned int C, MemType MEM, Centering CTR>
BoxData<T,C,MEM>& MBBoundaryRegister<T,C,MEM,CTR>::adjacent(
        MBIndex         a_localIndex,
        MBIndex         a_adjIndex)
{
    Box dirs = Box::Kernel(1);
    for (auto dir : dirs)
    {
        for (auto bound : bounds(a_localIndex, dir))
        {
            if (bound.adjIndex == a_adjIndex)
            {
                return *(bound.adjData);
            }
        }
    }
    MayDay<void>::Error("MBBoundaryRegister::adjacent | Error: no boundary between patches.");
}

// Get Rotation
template<typename T, unsigned int C, MemType MEM, Centering CTR>
CoordPermutation MBBoundaryRegister<T,C,MEM,CTR>::toLocal(
        MBIndex         a_index,
        Point           a_dir,
        unsigned int    a_adjBlock)
{
    return bound(a_index, a_dir, a_adjBlock).adjToLocal;
}

// Get Rotation
template<typename T, unsigned int C, MemType MEM, Centering CTR>
CoordPermutation MBBoundaryRegister<T,C,MEM,CTR>::toLocal(
        MBIndex         a_index,
        unsigned int    a_dir,
        Side::LoHiSide  a_side)
{
    return bound(a_index, a_dir, a_side).adjToLocal;
}

// Set all data to zero
template<typename T, unsigned int C, MemType MEM, Centering CTR>
void MBBoundaryRegister<T,C,MEM,CTR>::clearLocal()
{
    for (auto src : m_bounds)
    {
        for (auto bound : src)
        {
            bound.localData->setVal(0);
        }
    }
}

// Set all data to zero
template<typename T, unsigned int C, MemType MEM, Centering CTR>
void MBBoundaryRegister<T,C,MEM,CTR>::clearAdjacent()
{
    for (auto src : m_bounds)
    {
        for (auto bound : src)
        {
            bound.adjData->setVal(0);
        }
    }
}
// Set all data to zero
template<typename T, unsigned int C, MemType MEM, Centering CTR>
void MBBoundaryRegister<T,C,MEM,CTR>::clear()
{
    for (auto src : m_bounds)
    {
        for (auto& bound : src.second)
        {
            bound.localData->setVal(0);
            bound.adjData->setVal(0);
        }
    }
}
template<typename T, unsigned int C, MemType MEM, Centering CTR>
void MBBoundaryRegister<T,C,MEM,CTR>::
increment(Array<BoxData<T,C,MEM>,DIM>& a_fluxes,
        MBIndex a_dit,
        T a_weight)
{
    PR_TIMERS("MBBoundaryRegister::increment");
    auto layout = this -> layout();
    auto block = layout.block(a_dit);
    Box bx = layout[a_dit];
    for (auto mb : bounds(a_dit))
    {
        if (mb.localDir.codim() != 1) { continue; }
        int dir = mb.localFace().dir;
        auto side = mb.localFace().side;
        Box bxfaces = bx.extrude(Point::Basis(dir));
        Box bxLocal = bxfaces.face(dir,side) & a_fluxes[dir].box(); 
        BoxData<T,C,MEM> tmp = forall<T,C,MEM>
            ([ ] PROTO_LAMBDA(
                              Var<T,C,MEM>& a_corr,
                              const Var<T,C,MEM>& a_fluxdir,
                              T a_weight)
             {
                for (int comp = 0; comp < C; comp++)
                {
                    a_corr(comp) = a_fluxdir(comp)*a_weight;
                }
             },bxLocal,a_fluxes[dir],a_weight);
        if (!((mb.localData->box()).containsBox(bxLocal)))
        {
            cout << mb.localData -> box() << ", " << bxLocal << endl;
        }
        *mb.localData += tmp;
    } 
}
template<typename T, unsigned int C, MemType MEM, Centering CTR>
void MBBoundaryRegister<T,C,MEM,CTR>::
reflux(MBLevelBoxData<T,C,MEM,Centering::PR_CELL>& a_rhsTotal)
{
    PR_TIMERS("MBBoundaryRegister::reflux");
    auto layout = a_rhsTotal.layout();
    this->exchange();
    auto domain = layout.domain();
    for (auto dit : layout)
    {   
        auto block = layout.block(dit);
        auto nvec = layout.domain().getBlock(block).box().high();
        auto& a_rhsTotal_i = a_rhsTotal[dit];
        auto bx = layout[dit];
        for (auto mb : bounds(dit))
        {
            if (mb.localDir.codim() != 1) { continue; } //ignore non-face boundaries
            int dir = mb.localFace().dir;
            auto side = mb.localFace().side;
            T dx = 1.0/(nvec[dir]+1);
            auto R = mb.adjToLocal;
            auto Rinv = R.inverse();
            Point normal = Point::Basis(dir,side);
            Point adjNormal = Rinv(normal);
            auto adjblock = mb.adjIndex.local();          
            Point signflipPoint = Rinv(Point::Basis(dir));
            int signflip = signflipPoint.dot(Point::Ones());
            Array<int,C> signflipArray;
            signflipArray.fill(1);
            //#ifdef PR_CFD
            signflipArray.fill(signflip);
            //signflipArray[0] = 1;
            //signflipArray[4] = 1;
            //#endif
            PROTO_ASSERT(signflip*signflip == 1,"signflip from Coordinate transform != +-1");
            adjNormal *= -1;
            Point adjBasis = adjNormal*adjNormal;
            int dp = adjNormal.dot(adjBasis);
            Box bxfaces = bx.extrude(Point::Basis(dir));
            Box bxedge = bxfaces.face(dir,side);
            BoxData<double,C,HOST> adjCopy(bxedge);
            if (dp*(2*side-1) == 1)
            {
                mb.adjData->shift(Point::Basis(dir));
                mb.adjData -> copyTo(adjCopy);
                mb.adjData->shift(Point::Basis(dir,-1));
            }
            else
            {
                mb.adjData -> copyTo(adjCopy);
            }
            signflip = 1;
            BoxData<T,C,MEM> fluxCorr = forall<T,C,MEM>
              ([ ] PROTO_LAMBDA
               (Var<T,C,MEM>& a_corr,
                Var<T,C,MEM>& a_local,
                Var<T,C,MEM>& a_adj,
                Array<int,C> a_signfliparray)
               {
                 for (int comp = 0; comp < C; comp++)
                   { 
                     a_corr(comp) =
                     .5*(a_local(comp) - a_signfliparray[comp]*a_adj(comp));
                   }                       
               },*mb.localData,adjCopy,signflipArray);
            
            Stencil<T> faceStencil;
            if (side == Side::LoHiSide::Lo)
            {
                faceStencil = 1.0*Shift(Point::Zeros());
            }
            else if (side == Side::LoHiSide::Hi)
            { 
                faceStencil = (-1.0)*Shift(Point::Basis(dir));
            }          
            T wgt = -1.0/dx;
            BoxData<T,C,MEM> tmp = faceStencil(fluxCorr,wgt);
            a_rhsTotal_i += tmp;
        }
    }
}
// Exchange and fill adjacent data
/** Local data is maintained
*/
template<typename T, unsigned int C, MemType MEM, Centering CTR>
void MBBoundaryRegister<T,C,MEM,CTR>::exchange()
{
    PR_TIME("MBBoundaryRegister::exchange");
#ifdef PR_MPI
    {
        PR_TIME("MBBoundaryRegister::exchange (MPI barrier)");
        barrier();
    }
#endif
    m_copier.execute();
}

//=======================================================================================
// COPIER

template<typename T, unsigned int C, MemType MEM, Centering CTR>
MBBoundaryCopier<T,C,MEM,CTR>::MBBoundaryCopier(){}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void
MBBoundaryCopier<T,C,MEM,CTR>::buildMotionPlans(
        MBBoundaryCopierOp<T, C, MEM, CTR>& a_op)
{
    auto  layout = a_op.m_data->layout();
    auto  ghost = a_op.m_data->ghost();
    auto& graph = layout.domain().graph();
    Box dirs = Box::Kernel(1);
    int localRank = procID();

    for (auto iter : layout)
    {
        auto localBlock = layout.block(iter);
        for (auto dir : dirs)
        {
            auto bounds = a_op.m_data->bounds(iter, dir);
            for (auto bound : bounds)
            {
                PROTO_ASSERT(layout.procID(bound.localIndex) == localRank,
                        "MBBoundaryCopier::buildMotionPlans | Error: Data corruption.");
                int adjRank = layout.procID(bound.adjIndex);
                auto adjBlock = layout.block(bound.adjIndex);
                Point adjDir = graph.reverseDir(localBlock, adjBlock, dir);
                
                Box localBox_to = bound.localData->box();
                Box adjBox_to = layout.domain().convertBox(
                        bound.adjData->box(), localBlock, adjBlock);

                MBMotionItem toMotionItem(
                        bound.adjIndex, bound.localIndex, adjBox_to, localBox_to);
                toMotionItem.rotation = bound.adjToLocal;
                if (localRank == adjRank)
                {
                    this->m_localMotionPlan.push_back(toMotionItem);
                } else {
                    toMotionItem.procID = adjRank;
                    this->m_toMotionPlan.push_back(toMotionItem);
                  
                    Box localBox_from = localBox_to;
                    Box adjBox_from = adjBox_to;
                    
                    MBMotionItem fromMotionItem(
                            bound.localIndex, bound.adjIndex, localBox_from, adjBox_from);
                    fromMotionItem.rotation = bound.adjToLocal.inverse();
                    fromMotionItem.procID = adjRank;
                    this->m_fromMotionPlan.push_back(fromMotionItem);
                }
            }
        }
    }
    this->sort();
}

//=======================================================================================
// COPIER OP

template<typename T, unsigned int C, MemType MEM, Centering CTR>
MBBoundaryCopierOp<T,C,MEM,CTR>::MBBoundaryCopierOp(
        MBBoundaryRegister<T, C, MEM, CTR>& a_data)
{
    m_data = &a_data;
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
uint64_t
MBBoundaryCopierOp<T,C,MEM,CTR>::linearSize(
        const Box&          a_box,
        const MBIndex&      a_index) const
{
    return C * a_box.size() * sizeof(T);
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void
MBBoundaryCopierOp<T,C,MEM,CTR>::linearOut(
        void*               a_buffer,
        const MBMotionItem& a_info) const
{
    MemType mem = pointerMemType(a_buffer);
    PROTO_ASSERT(mem == MEM,
        "MBLevelExchangeCopierOp::linearOut | Error: Memtype mismatch. %i != %i", mem, MEM);
    const auto& srcIndex = a_info.fromIndex;
    const auto& dstIndex = a_info.toIndex;
    auto& patch = m_data->local(srcIndex, dstIndex); 
    CInterval cint(0, C - 1);
    patch.linearOut(a_buffer, a_info.fromRegion, cint);
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void
MBBoundaryCopierOp<T,C,MEM,CTR>::linearIn(
        void*               a_buffer,
        const MBMotionItem& a_info) const
{
    MemType mem = pointerMemType(a_buffer);
    PROTO_ASSERT(mem == MEM,
        "MBLevelExchangeCopierOp::linearOut | Error: Memtype mismatch. %i != %i", mem, MEM);
    const auto& srcIndex = a_info.fromIndex;
    const auto& dstIndex = a_info.toIndex;
    auto R = a_info.rotation;
    auto& patch = m_data->adjacent(dstIndex, srcIndex); 
    R.rotateBuffer<T,C>((T*)a_buffer, patch.data(), a_info.fromRegion, a_info.toRegion);
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void
MBBoundaryCopierOp<T,C,MEM,CTR>::localCopy(
        const MBMotionItem& a_info) const
{
    const auto &domain = a_info.fromRegion;
    const auto &srcIndex = a_info.fromIndex;
    const auto &range = a_info.toRegion;
    const auto &dstIndex = a_info.toIndex;
    
    bool linearizationTest = false;
    if (linearizationTest)
    {
        PR_TIMERS("MBBoundaryCopierOp::localCopy (linearization test)");
        size_t sizeSrc = linearSize(domain, srcIndex);
        size_t sizeDst = linearSize(range,  dstIndex);
        if (sizeSrc != sizeDst)
        {
            MayDay<void>::Abort(
                    "Linearization test failure. dst and src have different sizes");
        }
        void *srcBuffer = proto_malloc<MEM>(sizeSrc);
        void *dstBuffer = proto_malloc<MEM>(sizeDst);
        {
            linearOut(srcBuffer, a_info);
            proto_memcpy<MEM, MEM>(srcBuffer, dstBuffer, sizeSrc);
            linearIn(dstBuffer, a_info);
        }
        proto_free<MEM>(dstBuffer);
        proto_free<MEM>(srcBuffer);

    } else {
        PR_TIMERS("MBBoundaryCopierOp::localCopy");
        m_data->local(srcIndex, dstIndex).copyTo(
                m_data->adjacent(dstIndex, srcIndex),
                a_info.rotation);
    }
}







