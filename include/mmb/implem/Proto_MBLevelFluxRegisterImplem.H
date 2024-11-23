template<typename T, unsigned int C, MemType MEM>
const std::vector<Register<T,C,MEM>>& 
MBLevelFluxRegisterTester<T,C,MEM>::getCoarseRegistersAtIndex(MBIndex a_index) const
{
    return m_register.m_crseRegisters[a_index];
}

template<typename T, unsigned int C, MemType MEM>
const std::vector<Register<T,C,MEM>>& 
MBLevelFluxRegisterTester<T,C,MEM>::getFineRegistersAtIndex(MBIndex a_index) const
{
    return m_register.m_fineRegisters[a_index];
}


template<typename T, unsigned int C, MemType MEM>
MBFluxRegisterCopierOp<T,C,MEM>::MBFluxRegisterCopierOp(MBLevelFluxRegister<T, C, MEM>& a_register)
{
    m_register = &a_register;
}

template<typename T, unsigned int C, MemType MEM>
int MBFluxRegisterCopierOp<T,C,MEM>::linearSize(
                                   const Box& a_bx,
                                   const MBIndex& a_index) const
{
    return C*a_bx.size()*sizeof(T);
}
template<typename T, unsigned int C, MemType MEM>
void MBFluxRegisterCopierOp<T,C,MEM>::linearOut(void* a_buf, const MBMotionItem& a_info) const
{
    // linearOut copies from source data holders into buffers.
    const auto& index = a_info.fromIndex;
    const auto& domain = a_info.fromRegion;
    auto& data = m_register->sourceData(domain, index);
    CInterval cint(0,C-1);
    data.linearOut(a_buf, domain, cint);
}

template<typename T, unsigned int C, MemType MEM>
void MBFluxRegisterCopierOp<T,C,MEM>::linearIn(void* a_buf, const MBMotionItem& a_info) const
{
    // linearIn copies from buffers into destination data holders.
    const auto& index = a_info.toIndex;
    const auto& range = a_info.toRegion;
    //BoxData<T,C,MEM>& data = m_register->destData(a_bx,a_index);
    auto& data = m_register->destData(range, index);
    CInterval cint(0,C-1);
    data.linearIn(a_buf, range, cint);
}

template<typename T, unsigned int C, MemType MEM>
void MBFluxRegisterCopierOp<T,C,MEM>::localCopy(const MBMotionItem& a_info) const
{
#define noLinearizationTest 1
    const auto& srcIndex =  a_info.fromIndex;
    const auto& domain =    a_info.fromRegion;
    const auto& dstIndex =  a_info.toIndex;
    const auto& range =     a_info.toRegion;
    
#if noLinearizationTest

    PR_TIMERS("MBFluxRegisterCopier_localCopy_1");
    auto& src = m_register->sourceData(domain, srcIndex);
    auto& dst = m_register->destData(  range,  dstIndex);
    Point shift = range.low() - domain.low();
    src.copyTo(dst, domain, shift);
#else
    // code for testing linearIn / linearOut while running serial.
    PR_TIMERS("MBFluxRegisterCopier_localCopy_2");
    int srcSize = linearSize(domain, srcIndex);
    int dstSize = linearSize(range,  dstIndex);

    if (srcSize != dstSize)
    {
        MayDay<void>::Abort("LinearizationTest failure: dest and source have different sizes");
    }

    std::vector<char> buffer(srcSize);
    {
        void* b = (void*)&(buffer[0]);
        linearOut(b, a_info);
        linearIn (b, a_info);
    }
#endif
}

template<typename T, unsigned int C, MemType MEM>
void MBFluxRegisterCopier<T,C,MEM>::buildMotionPlans (MBFluxRegisterCopierOp<T, C, MEM >& a_op)
{
    auto & lfr = *(a_op.m_register);
    auto & fineLayout = lfr.m_fineLayout;
    auto & crseLayout = lfr.m_crseLayout;
    auto & crseRegisters = lfr.m_crseRegisters;
    auto & crseIncrement = lfr.m_crseIncrement;
    auto & fineRegisters = lfr.m_fineRegisters;
    auto & fineIndices = lfr.m_fineIndices;

    for (auto citer : crseLayout)
    {
        auto block = crseLayout.block(citer);
        auto & refRatio = lfr.m_refRatios[block];
        auto& localRegisters = crseRegisters[citer];
        for (int k = 0; k < localRegisters.size(); k++)
        {
            // Find the index of the fineRegister source.
            // The index is associated with the fine patch contributing the averaged flux
            auto& reg     = localRegisters[k];
            Box   destBox = reg.m_data->box();
            auto  side    = reg.m_side;
            int   dir     = reg.m_dir;
            Point finePatchPoint = 
                destBox.low()*refRatio/fineLayout.boxSize(block) + Point::Basis(dir,side);
            auto fineIndex = fineLayout.find(finePatchPoint, block);
            
            PROTO_ASSERT(fineIndex != *fineLayout.end(),
                "MBFluxRegisterCopier::buildMotionPlans | Error:\
                Could not find fine register index.");
            
            // Compute fine register box (nontrivial b/c of possibility of periodic BCs. 
            Box cfBox  = fineLayout.box(fineIndex).coarsen(refRatio);
            Box srcBox = cfBox.adjacent(dir, flip(side), 1);
            
            // Boxes for the fine registers are defined relative to the fine patch
            // coordinates, even if the register Boxes are periodic images outside the
            // ProblemDomain Box.
            MBMotionItem toMotionItem(fineIndex, citer, srcBox, destBox);
            auto fineID = fineLayout.procID(fineIndex);
            if (fineID == procID())
            { 
                this->m_localMotionPlan.push_back(toMotionItem);
            }
            else
            {
                toMotionItem.procID = fineID;
                this->m_toMotionPlan.push_back(toMotionItem);
            }
        } // for each local coarse register
    } // end for coarse patch

    // Iterate over fine grids that I own to obtain entries in m_fromMotionPlan.           
    int destProcID = procID();
    for (auto fiter : fineLayout.begin())
    {
        auto block = fineLayout.block(fiter);
        for (int dir = 0; dir < DIM; dir++)
        {
            SideIterator siter;
            for (siter.begin(); siter.ok(); ++siter)
            {
                int k = fineIndices[fiter][dir][*siter];
                if (k >= 0) //otherwise there is no fine register for this (dir, side)
                {
                  Box srcBox = fineRegisters[fiter][k].m_data->box();
                  Point crsePatchPoint = srcBox.low()/crseLayout.boxSize(block);
                    
                    auto crseIndex  = crseLayout.find(crsePatchPoint, block);
                    int  crseProcID = crseLayout.procID(crseIndex);
                    //Automatically defines dstBox to be the correct periodic image. NO IT DOESN'T.
                    //Box dstBox = crseLayout.domain() & srcBox;
                    //dstBox = crseLayout.domain().image(dstBox); 
                    
                    // TODO: consider periodic
                    Box dstBox = srcBox;
                    if (crseProcID != procID())
                    {
                        MBMotionItem fromMotionItem(fiter, crseIndex, srcBox, dstBox);
                        fromMotionItem.procID = crseProcID;
                        this->m_fromMotionPlan.push_back(fromMotionItem);
                    }
                }
            }
        }
    }     
    this->sort();
}

template<typename T, unsigned int C, MemType MEM>
MBLevelFluxRegister<T,C,MEM>::MBLevelFluxRegister(
        const MBDisjointBoxLayout& a_crseLayout,
        const MBDisjointBoxLayout& a_fineLayout,
        const std::vector<Point>&   a_refRatios,
        std::vector<Array<T, DIM>>& a_dxCoarse)
{
    define(a_crseLayout,a_fineLayout,a_refRatios,a_dxCoarse);
}

template<typename T, unsigned int C, MemType MEM>
MBLevelFluxRegister<T,C,MEM>::MBLevelFluxRegister(
        const MBDisjointBoxLayout& a_crseLayout,
        const MBDisjointBoxLayout& a_fineLayout,
        Point&          a_refRatio,
        Array<T, DIM>&  a_dxCoarse)
{
    define(a_crseLayout, a_fineLayout, a_refRatio, a_dxCoarse);
}

template<typename T, unsigned int C, MemType MEM>
void MBLevelFluxRegister<T,C,MEM>::define (
        const MBDisjointBoxLayout& a_crseLayout,
        const MBDisjointBoxLayout& a_fineLayout,
        Point&          a_refRatio,
        Array<T, DIM>&  a_dxCoarse)
{
    int numBlocks = a_crseLayout.numBlocks();
    std::vector<Point> refRatios(numBlocks, a_refRatio);
    std::vector<Array<T, DIM>> dx(numBlocks, a_dxCoarse);
    define(a_crseLayout, a_fineLayout, refRatios, dx);
}


template<typename T, unsigned int C, MemType MEM>
void MBLevelFluxRegister<T,C,MEM>::define (
        const MBDisjointBoxLayout& a_crseLayout,
        const MBDisjointBoxLayout& a_fineLayout,
        const std::vector<Point>&   a_refRatios,
        std::vector<Array<T, DIM>>& a_dxCoarse)
{

    validateInputs(a_crseLayout, a_fineLayout, a_refRatios, a_dxCoarse);

    m_refRatios = a_refRatios;
    m_dxCrse = a_dxCoarse;
    m_fineLayout = a_fineLayout;
    m_crseLayout = a_crseLayout;

    clearInternalData();
    defineCoarse();
    defineFine();
    defineAverageStencils();
    
    MBFluxRegisterCopierOp<T,C,MEM> op(*this);
    m_copier.define(op);

    initCoarseRegisters(0);
    initFineRegisters(0);
    initIncrement(0);
}

template<typename T, unsigned int C, MemType MEM>
void MBLevelFluxRegister<T,C,MEM>::validateInputs (
        const MBDisjointBoxLayout& a_crseLayout,
        const MBDisjointBoxLayout& a_fineLayout,
        const std::vector<Point>&   a_refRatios,
        std::vector<Array<T, DIM>>& a_dxCoarse) const
{
    int numBlocks = a_crseLayout.numBlocks();
    PROTO_ASSERT(a_fineLayout.numBlocks() == numBlocks,
        "MBLevelFluxRegister::validateInputs | Error: layouts have different block counts.");
    PROTO_ASSERT(a_refRatios.size() == numBlocks,
        "MBLevelFluxRegister::validateInputs | Error: wrong number of refinement ratios.");
    PROTO_ASSERT(a_dxCoarse.size() == numBlocks,
        "MBLevelFluxRegister::validateInputs | Error: wrong number of grid spacings.");

    for (int bi = 0; bi < numBlocks; bi++)
    {
        Point ri = a_refRatios[bi];
        Point patchRefRatio = 
            a_crseLayout.boxSize(bi)*ri / a_fineLayout.boxSize(bi);
        PROTO_ASSERT(patchRefRatio*a_fineLayout.boxSize(bi) == a_crseLayout.boxSize(bi)*ri,
            "MBLevelFluxRegister::validateInputs | Error: coarse region is not covered by fine region.");
        
        Box crseDomainBox = a_crseLayout.domain().getBlock(bi).box();
        Box fineDomainBox = a_fineLayout.domain().getBlock(bi).box();
        PROTO_ASSERT(fineDomainBox.containsBox(crseDomainBox),
            "MBLevelFluxRegister::validateInputs | Error: Coarse and fine layouts appear to be reversed");
    }

    // Check to see whether coarsened fine patch is more than one coarse cell wide in all directions that are being refined. 
    for (int bi = 0; bi < numBlocks; bi++)
    {
        for (int dir = 0; dir < DIM; dir++)
        {
            Point ri = a_refRatios[bi];
            if (ri[dir] > 1)
            {
                PROTO_ASSERT(a_fineLayout.boxSize(bi)[dir]/ri[dir] > 1,
                    "MBLevelFluxRegister::validateInputs | Error:\
                    fine patch cell size is smaller than the refinement ratio.");
            }
        }
    }
}

template<typename T, unsigned int C, MemType MEM>
void MBLevelFluxRegister<T,C,MEM>::clearInternalData ()
{
    m_crseRegisters.clear();
    m_crseIncrement.clear();
    m_crseIndices.clear();
    m_fineRegisters.clear();
    m_fineIndices.clear();

    int numCrsePatches = m_crseLayout.localSize();

    m_crseRegisters.resize(numCrsePatches);
    m_crseIncrement.resize(numCrsePatches);
    m_crseIndices.resize(numCrsePatches);
}


template<typename T, unsigned int C, MemType MEM>
void MBLevelFluxRegister<T,C,MEM>::defineCoarse ()
{
    auto& crsePatchDomain = m_crseLayout.partition().domain();
    auto& finePatchDomain = m_fineLayout.partition().domain();

    for (auto citer : m_crseLayout)
    {
        m_crseRegisters[citer].clear();
        m_crseIncrement[citer].clear();
        m_crseIndices[citer].clear();

        BlockIndex block = m_crseLayout.block(citer);
        Point refRatio = m_refRatios[block];
        Point fineBoxSize = m_fineLayout.boxSize(block);
        Box crseBox = m_crseLayout[citer];

        auto coveringFinePatches = crseBox.refine(refRatio).coarsen(fineBoxSize);
        PROTO_ASSERT(finePatchDomain.getBlock(block).box().containsBox(coveringFinePatches),
            "MBLevelFluxRegister::define | Error: Data corruption.");
        for (auto finePatch : coveringFinePatches)
        {
            addCoarseRegisters(finePatch, citer);
        }
    }
/*
    for (auto citer : m_crseLayout)
    {
        BlockIndex block = m_crseLayout.block(citer);
        Point refRatio = m_refRatios[block];
        Point fineBoxSize = m_fineLayout.boxSize(block);
        Point crseBoxSize = m_crseLayout.boxSize(block);
        Box crseBox = m_crseLayout[citer];
        m_crseRegisters[citer].clear();
        m_crseIncrement[citer].clear();
        m_crseIndices[citer].clear();

        auto coveringFinePatches = crseBox.refine(refRatio).coarsen(fineBoxSize);
        PROTO_ASSERT(finePatchDomain.getBlock(block).box().containsBox(coveringFinePatches),
            "MBLevelFluxRegister::define | Error: Data corruption.");
        int k = 0;
        for (auto finePatch : coveringFinePatches)
        {
            // check if fine patch is ABSENT from fine layout
            auto fineIndex = m_fineLayout.find(finePatch, block);
            if (fineIndex != *m_fineLayout.end()) { continue; }
            // check adjacent fine patch that is PRESENT in fine layout
            Box cfBox = Box(finePatch, finePatch).refine(fineBoxSize).coarsen(refRatio);
            for (auto dir : Point::DirectionsOfCodim(1))
            {
                MBIndex fineIndexAdj = findAdjacentFinePatchIndex(fineIndex, finePatch, dir);
                if (fineIndexAdj == *m_fineLayout.end()) { continue; } //adj patch is absent
                // found a valid cf boundary in block interior
                Box destBox = cfBox.face(dir);

                int thisKey = computeTileIndex(destBox, citer);
                shared_ptr<BoxData<T,C,MEM> > crseRegisterDataPtr(new BoxData<T,C,MEM>(destBox));
                shared_ptr<BoxData<T,C,MEM> > crseIncrementDataPtr(new BoxData<T,C,MEM>(destBox));
                Register<T,C,MEM> reg(crseRegisterDataPtr, dir);
                m_crseRegisters[citer].push_back(reg);
                m_crseIncrement[citer].push_back(crseIncrementDataPtr);
                std::cout << "k: " << k << " | size: " << m_crseRegisters[citer].size() << std::endl;
                m_crseIndices  [citer][thisKey] = k;
                k++;
            } // end dir loop
        } // end fine patch loop
    } // end coarse loop
    */
}

template<typename T, unsigned int C, MemType MEM>
void MBLevelFluxRegister<T,C,MEM>::addCoarseRegisters(
    Point a_finePatch,
    MBIndex a_coarseIndex)
{
    BlockIndex block = m_crseLayout.block(a_coarseIndex);
    auto fineIndex = m_fineLayout.find(a_finePatch, block);

    // this region must be unrefined for a coarse-fine boundary to exist
    // we assume THIS fine patch is unrefined and search for adjacent refined patches
    if (fineIndex != *m_fineLayout.end()) { return; }

    for (auto dir : Point::DirectionsOfCodim(1))
    {
        MBIndex fineIndexAdj = findAdjacentFinePatchIndex(fineIndex, a_finePatch, dir);
        addCoarseRegisterFromAdjacent(a_finePatch, a_coarseIndex, dir, fineIndexAdj);
    }
}

template<typename T, unsigned int C, MemType MEM>
MBIndex MBLevelFluxRegister<T,C,MEM>::findAdjacentFinePatchIndex(
    MBIndex a_fineIndex,
    Point a_finePatch,
    Point a_dir) const
{
    auto& finePatchDomain = m_fineLayout.partition().domain();
    BlockIndex block = m_fineLayout.block(a_fineIndex);
    Point adjacentPatch = a_finePatch + a_dir;
    MBIndex adjacentIndex = *m_fineLayout.end();
    if (finePatchDomain.inInterior(adjacentPatch, block))
    {
        // INTERIOR case (normal case)
        adjacentIndex = m_fineLayout.find(adjacentPatch, block);
    } else if (finePatchDomain.inBlockBoundary(adjacentPatch, block))
    {
        // BLOCK BOUNDARY case
        BlockIndex adjacentBlock = finePatchDomain.graph().adjacent(block, a_dir);
        Point adjacentPatchLocalCoords = 
            finePatchDomain.convert(adjacentPatch, block, adjacentBlock);
        adjacentIndex = m_fineLayout.find(adjacentPatchLocalCoords, adjacentBlock);
    }
    else if (finePatchDomain.inDomainBoundary(adjacentPatch, block))
    {
        // EXTERIOR case (no CF boundary)
    } else {
        // INVALID case
        PROTO_ASSERT(false, "MBLevelFluxRegister::define | Error: Data corruption.");
    }  
    return adjacentIndex;
}

template<typename T, unsigned int C, MemType MEM>
void MBLevelFluxRegister<T,C,MEM>::addCoarseRegisterFromAdjacent(
    Point a_finePatch,
    MBIndex a_coarseIndex,
    Point a_dir,
    MBIndex a_fineIndexAdj)
{
    // we are specifically searching for refined regions since the
    // patch we are searching from (a_finePatch) is unrefined
    if (a_fineIndexAdj == *m_fineLayout.end()) { return; }
    
    BlockIndex block = m_crseLayout.block(a_coarseIndex);
    Point refRatio = m_refRatios[block];
    Point fineBoxSize = m_fineLayout.boxSize(block);
    int numCoarseRegisters = m_crseRegisters[a_coarseIndex].size();

    Box cfBox = Box(a_finePatch, a_finePatch).refine(fineBoxSize).coarsen(refRatio);
    Box destBox = cfBox.face(a_dir);

    shared_ptr<BoxData<T, C, MEM>> crseRegisterDataPtr(new BoxData<T, C, MEM>(destBox));
    Register<T, C, MEM> reg(crseRegisterDataPtr, a_dir);
    m_crseRegisters[a_coarseIndex].push_back(reg);

    shared_ptr<BoxData<T, C, MEM>> crseIncrementDataPtr(new BoxData<T, C, MEM>(destBox));
    m_crseIncrement[a_coarseIndex].push_back(crseIncrementDataPtr);

    int tileID = computeTileIndex(destBox, a_coarseIndex);
    m_crseIndices[a_coarseIndex][tileID] = numCoarseRegisters;
}

// a_bx is the coarse box updated by refluxing
template<typename T, unsigned int C, MemType MEM>
int MBLevelFluxRegister<T,C,MEM>::computeTileIndex(const Box& a_bx, const MBIndex& a_di)
{
    // Find which direction is normal to the register. We use the fact that the register
    // is of size 1 only in the normal direction.
   
    Box crseBox = m_crseLayout[a_di];
    Point origin = crseBox.low();
    Box tileBox = a_bx.shift(-origin);
    Point sizes = crseBox.high() - crseBox.low() + Point::Ones();

    int k = tileBox.low()[0];
    int p = 1;
    for (int ii = 1; ii < DIM; ii++)
    {
        p *= sizes[ii-1];
        k += p*tileBox.low()[ii];
    }
    p *= sizes[DIM-1];
    k += p*tileBox.high()[0];
    for (int ii = 1; ii < DIM; ii++)
    {
        p *= sizes[ii-1];
        k += p*tileBox.high()[ii];
    }
    
    return k;
}

template<typename T, unsigned int C, MemType MEM>
void MBLevelFluxRegister<T,C,MEM>::defineFine ()
{
    auto& crsePatchDomain = m_crseLayout.partition().domain();
    auto& finePatchDomain = m_fineLayout.partition().domain();
    
    // Set up fine data holders.
    int numFinePatches = m_fineLayout.localSize();
    m_fineRegisters.resize(numFinePatches);
    m_fineIndices.resize(numFinePatches);
    for (auto fiter : m_fineLayout)
    {
        BlockIndex block = m_fineLayout.block(fiter);
        Point finePoint = m_fineLayout.point(fiter);
        Box fineBox = m_fineLayout[fiter];
        Box cfBox = fineBox.coarsen(m_refRatios[block]);
        int k = 0;
        for (auto dd = 0; dd < DIM; dd++)
        {
            SideIterator siter;
            for (siter.begin(); siter.ok(); ++siter)
            {
                Point dir = Point::Basis(dd,*siter);
                Point finePointAdj = finePoint + dir;
                auto fineIndexAdj = *m_fineLayout.end();
                if (finePatchDomain.inInterior(finePointAdj, block))
                {
                    // INTERIOR case (normal case)
                    fineIndexAdj = m_fineLayout.find(finePointAdj, block);
                } else if (finePatchDomain.inBlockBoundary(finePointAdj, block))
                {
                    // BLOCK case
                    BlockIndex blockAdj = finePatchDomain.graph().adjacent(block, dir);
                    Point finePointAdjLocal = 
                        finePatchDomain.convert(finePointAdj, block, blockAdj);
                    fineIndexAdj = m_fineLayout.find(finePointAdjLocal, blockAdj);
                }
                else if (finePatchDomain.inDomainBoundary(finePointAdj, block))
                {
                    // EXTERIOR case (no CF boundary)
                    // TODO: do I need this? m_fineIndices[fiter][d][(int)(*siter)] = -1;
                    continue;
                } else {
                    // INVALID case
                    PROTO_ASSERT(false,
                            "MBLevelFluxRegister::define | Error: Data corruption.");
                }
                if (fineIndexAdj == *m_fineLayout.end())
                {
                    //found a cf boundary between finePoint and finePointAdj
                    Box srcBox = cfBox.adjacent(dir);

                    shared_ptr<BoxData<T,C,MEM> > temp(new BoxData<T,C,MEM>(srcBox));
                    Register<T,C,MEM> reg(temp, dd, *siter);
                    m_fineRegisters[fiter].push_back(reg);
                    m_fineIndices[fiter][dd][(int)(*siter)] = k;
                    k++;
                } else {
                    m_fineIndices[fiter][dd][(int)(*siter)] = -1;
                }
            }
        }
    }
}


template<typename T, unsigned int C, MemType MEM>
void MBLevelFluxRegister<T,C,MEM>::defineAverageStencils ()
{
    m_avgFlux.resize(m_crseLayout.numBlocks());
    for (BlockIndex bi = 0; bi < m_crseLayout.numBlocks(); bi++)
    {
        for (int d = 0; d < DIM; d++)
        {
            SideIterator siter;
            for (siter.begin(); siter.ok(); ++siter)
            {
                // Side::Lo is not a typo -CLG
                auto avg = Stencil<T>::AvgDownFace(d, Side::Lo, m_refRatios[bi]);
                if (*siter == Side::Hi)
                {
                    avg = avg * (-1);
                }
                else if (*siter == Side::Lo)
                {
                    // Side::Lo is not a typo here either -CLG
                    avg.destShift() = Point::Basis(d, Side::Lo);
                }
                m_avgFlux[bi][d][*siter] = avg;
            }
        }
    }
}

template<typename T, unsigned int C, MemType MEM>
void MBLevelFluxRegister<T,C,MEM>::incrementCoarseRegister (
        const BoxData<T,C,MEM>& a_flux,
        const MBIndex& a_crseIndex,
        const T& a_weight,
        unsigned int a_dir)
{
    auto& crseRegisters = m_crseRegisters[a_crseIndex];
    T weight = a_weight / m_dxCrse[a_dir];
    Stencil<T> scaleHigh = weight*Shift::Basis(a_dir, +1);
    Stencil<T> scaleLow = -weight*Shift::Zeros();
    for (int k = 0; k < crseRegisters.size(); k++)
    {   
        auto& reg = crseRegisters[k];
        if (reg.m_dir == a_dir)
        {
            auto& data = reg.m_data; 
            auto  side = reg.m_side;
            if (side == Side::Lo)
            {
                (*data) += scaleLow(a_flux);
            }
            else if (side == Side::Hi)
            {
                (*data) += scaleHigh(a_flux);
            }
        }
    }
}

template<typename T, unsigned int C, MemType MEM>
void MBLevelFluxRegister<T,C,MEM>::incrementFineRegister (
        const BoxData<T, C , MEM>& a_flux,
        const MBIndex& a_fineIndex,
        const T& a_weight,
        unsigned int a_dir)
{
    auto & fineRegister = m_fineRegisters[a_fineIndex];
    T weight = a_weight / m_dxCrse[a_dir];
    for (int k = 0; k < fineRegister.size(); k++)
    {
        auto& reg = fineRegister[k];
        if (reg.m_dir == a_dir)
        {
            auto& data = reg.m_data;
            auto  side = reg.m_side;
            (*data) += m_avgFlux[a_dir][side](a_flux, weight);
        }
    }
}

template<typename T, unsigned int C, MemType MEM>
void MBLevelFluxRegister<T,C,MEM>::applyRefluxCorrection(
        MBLevelBoxData<T,C,MEM>& a_coarseData,
        const T& a_weight)
{
    PROTO_ASSERT(m_crseLayout.compatible(a_coarseData.layout()),
            "MBLevelFluxRegister::reflux | Error: \
            Input data has incompatible layout.");
   
    // the copier moves fine register data into the coarse increment. 
    m_copier.execute();

    // We compute the difference between the coarse and fine registers, and increment a_coarseData.
    Stencil<T> increment = (-1.0)*Shift::Zeros();
    for (auto iter = a_coarseData.begin(); iter.ok(); ++iter)
    {
        int numRegisters = m_crseIncrement[*iter].size();
        auto& dest = m_crseIncrement[*iter];
        auto& creg = m_crseRegisters[*iter];
        // subtract off the coarse fluxes
        for (int k = 0; k < numRegisters; k++)
        {
            *(dest[k]) += increment(*(creg[k].m_data));
            *(dest[k]) *= a_weight;
        }
        // add increment to input data
        for (int k = 0; k < numRegisters; k++)
        {
            a_coarseData[*iter] += *(dest[k]);
        }
    }          
}

template <typename T, unsigned int C, MemType MEM>
void MBLevelFluxRegister<T, C, MEM>::initCoarseRegisters(T a_val)
{
    for (auto citer = m_crseLayout.begin(); citer.ok(); ++citer)
    {
        for (int k = 0; k < m_crseIndices[*citer].size(); k++)
        {
            m_crseRegisters[*citer][k].m_data->setVal(a_val);
        }
    }
}

template <typename T, unsigned int C, MemType MEM>
void MBLevelFluxRegister<T, C, MEM>::initFineRegisters(T a_val)
{
    for (auto fiter = m_fineLayout.begin(); fiter.ok(); ++fiter)
    {
        for (int k = 0; k < m_fineRegisters[*fiter].size(); k++)
        {
            m_fineRegisters[*fiter][k].m_data->setToZero();
        }
    }
}

template <typename T, unsigned int C, MemType MEM>
void MBLevelFluxRegister<T, C, MEM>::initIncrement(T a_val)
{
    for (auto citer = m_crseLayout.begin(); citer.ok(); ++citer)
    {
        for (int k = 0; k < m_crseIndices[*citer].size(); k++)
        {
            m_crseIncrement[*citer][k]->setToZero();
        }
    }
}



template<typename T, unsigned int C, MemType MEM>
BoxData<T,C,MEM>& MBLevelFluxRegister<T,C,MEM>::sourceData(const Box& a_bx, const MBIndex& a_fineIndex)
{
    auto block = m_fineLayout.block(a_fineIndex);
    std::shared_ptr<BoxData<T,C,MEM>> retval;
    Box fineBox = m_fineLayout.box(a_fineIndex);
    Box cfBox = fineBox.coarsen(m_refRatios[block]);
    SideIterator siter;
    for (siter.begin(); siter.ok(); ++siter)
    {
        for (int dir = 0; dir < DIM; dir++)
        {
            if (a_bx == cfBox.adjacent(dir, *siter, 1))
            {
                int k = m_fineIndices[a_fineIndex][dir][(int)(*siter)];
                retval = m_fineRegisters[a_fineIndex][k].m_data;
            }
        }
    }
    if (!retval) {
        std::cout << "Invalid Box: " << a_bx << " or invalid index: " << (int)a_fineIndex << std::endl;
        PROTO_ASSERT(false,
            "MBLevelFluxRegister::sourceData | Error: Could not locate data.");
    }
    return *retval;
}

template<typename T, unsigned int C, MemType MEM>
BoxData<T,C,MEM>& MBLevelFluxRegister<T,C,MEM>::destData(const Box& a_bx, const MBIndex& a_crseIndex)
{
    int crseKey       = computeTileIndex(a_bx, a_crseIndex);
    int mapIndex      = m_crseIndices[a_crseIndex][crseKey];
    return *(m_crseIncrement[a_crseIndex][mapIndex]);
}

template<typename T, unsigned int C, MemType MEM>
void MBLevelFluxRegister<T,C,MEM>::print() const
{
    using Proto::pr_out;
    pr_out() << "Level Flux Register: " << std::endl;
    pr_out() << "Coarse Registers: " << std::endl;
    for (int ii = 0; ii < m_crseRegisters.size(); ii++)
    {
        for (int jj = 0; jj < m_crseRegisters[ii].size(); jj++)
        {
            auto item = m_crseRegisters[ii][jj];
            pr_out() << "\t box: " << item.m_data->box();
            pr_out() << " | dir: " << item.m_dir;
            pr_out() << " | side: " << item.m_side << std::endl;
        }
    }
    pr_out() << "Fine Registers: " << std::endl;
    for (int ii = 0; ii < m_fineRegisters.size(); ii++)
    {
        for (int jj = 0; jj < m_fineRegisters[ii].size(); jj++)
        {
            auto item = m_fineRegisters[ii][jj];
            pr_out() << "\t box: " << item.m_data->box();
            pr_out() << " | dir: " << item.m_dir;
            pr_out() << " | side: " << item.m_side << std::endl;
        }
    }
}
