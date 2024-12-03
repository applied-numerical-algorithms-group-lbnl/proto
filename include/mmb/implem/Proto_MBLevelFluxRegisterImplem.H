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

template <typename T, unsigned int C, MemType MEM>
int MBFluxRegisterCopierOp<T, C, MEM>::linearSize(
    const Box &a_box,
    const MBIndex &a_index) const
{
    return C*a_box.size()*sizeof(T);
}
template<typename T, unsigned int C, MemType MEM>
void MBFluxRegisterCopierOp<T,C,MEM>::linearOut(void* a_buffer, const MBMotionItem& a_motionItem) const
{
    // linearOut copies from source data holders into buffers.
    const auto& index = a_motionItem.fromIndex;
    const auto& domain = a_motionItem.fromRegion;
    auto& data = m_register->sourceData(domain, index);
    CInterval cint(0,C-1);
    data.linearOut(a_buffer, domain, cint);
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
    auto & fluxRegister = *(a_op.m_register);

    addCoarseMotionItems(fluxRegister);
    addFineMotionItems(fluxRegister);
    this->sort();
}

template<typename T, unsigned int C, MemType MEM>
void MBFluxRegisterCopier<T,C,MEM>::addCoarseMotionItems (MBLevelFluxRegister<T,C,MEM>& a_fluxRegister)
{
    auto& coarseLayout = a_fluxRegister.m_crseLayout;
    auto& fineLayout = a_fluxRegister.m_fineLayout;
    for (auto iter : coarseLayout)
    {
        auto block = coarseLayout.block(iter);
        auto & refRatio = a_fluxRegister.m_refRatios[block];
        auto& localRegisters = a_fluxRegister.m_crseRegisters[iter];
        Point fineBoxSize = fineLayout.boxSize(block);

        for (auto& ri : localRegisters)
        {
            Box coarseBox = ri.data().box();
            Box fineBox = coarseBox.refine(refRatio);
            PatchID finePatch = fineBox.low() / fineBoxSize;
            MBIndex adjFineIndex = fineLayout.findAdjacentIndex(finePatch, block, ri.dir());
            PROTO_ASSERT(adjFineIndex != *fineLayout.end(),
                "MBFluxRegisterCopier::addCoarseMotionItems | Error: Data corruption");
            
            MBMotionItem toMotionItem(adjFineIndex, iter, coarseBox, coarseBox);
            auto fineID = fineLayout.procID(adjFineIndex);
            if (fineID == procID())
            { 
                this->m_localMotionPlan.push_back(toMotionItem);
            }
            else
            {
                toMotionItem.procID = fineID;
                this->m_toMotionPlan.push_back(toMotionItem);
            }
        }
    }
}

template<typename T, unsigned int C, MemType MEM>
void MBFluxRegisterCopier<T,C,MEM>::addFineMotionItems (MBLevelFluxRegister<T,C,MEM>& a_fluxRegister)
{
    auto& coarseLayout = a_fluxRegister.m_crseLayout;
    auto& fineLayout = a_fluxRegister.m_fineLayout;
    for (auto fineIndex : fineLayout.begin())
    {
        auto block = fineLayout.block(fineIndex);
        for (Point dir : Point::DirectionsOfCodim(1))
        {
            int coord = dir.firstNonZeroIndex();
            int side = (dir[coord] < 0) ? 0 : 1;
            int registerIndex = a_fluxRegister.m_fineIndices[fineIndex][coord][side];

            // register index of -1 means that there is no register associated with dir
            if (registerIndex < 0) { continue; }

            Register<T,C,MEM>& fineRegister = a_fluxRegister.m_fineRegisters[fineIndex][registerIndex];

            Box registerBox = fineRegister.data().box();

            // If this code needs to deal with periodicity, this will no longer be trivial
            // because one of these boxes may need to be shifted across the periodic boundary
            Box srcBox = registerBox; 
            Box dstBox = registerBox;

            MBIndex adjCoarseIndex = getCoarseIndex(a_fluxRegister, fineIndex, dir);
            int adjCoarseProcID = coarseLayout.procID(adjCoarseIndex);
            if (adjCoarseProcID != procID())
            {
                MBMotionItem fromMotionItem(fineIndex, adjCoarseIndex, srcBox, dstBox);
                fromMotionItem.procID = adjCoarseProcID;
                this->m_fromMotionPlan.push_back(fromMotionItem);
            }
        }
    }     
}

template<typename T, unsigned int C, MemType MEM>
MBIndex MBFluxRegisterCopier<T,C,MEM>::getCoarseIndex(
    MBLevelFluxRegister<T,C,MEM>& a_fluxRegister,
    MBIndex a_fineIndex,
    Point a_dir)
{
    auto& fineLayout = a_fluxRegister.m_fineLayout;
    auto& coarseLayout = a_fluxRegister.m_crseLayout;
    
    BlockIndex block = fineLayout.block(a_fineIndex);
    Point refRatio = a_fluxRegister.m_refRatios[block];

    PatchID finePatch = fineLayout.point(a_fineIndex);
    PatchID adjFinePatch = finePatch + a_dir;
    PatchID coarsePatch = finePatch / refRatio;
    PatchID adjCoarsePatch = adjFinePatch / refRatio;
    MBIndex coarseIndex = coarseLayout.find(coarsePatch, block);
    if (coarsePatch == adjCoarsePatch)
    {
        return coarseIndex;
    } else {
        return coarseLayout.findAdjacentIndex(coarsePatch, block, a_dir);
    }
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
}

template<typename T, unsigned int C, MemType MEM>
void MBLevelFluxRegister<T,C,MEM>::addCoarseRegisters(
    PatchID a_finePatch,
    MBIndex a_coarseIndex)
{
    BlockIndex block = m_crseLayout.block(a_coarseIndex);
    auto fineIndex = m_fineLayout.find(a_finePatch, block);

    // this region must be unrefined for a coarse-fine boundary to exist
    // we assume THIS fine patch is unrefined and search for adjacent patches that are refined
    if (fineIndex != *m_fineLayout.end()) { return; }

    for (auto dir : Point::DirectionsOfCodim(1))
    {
        MBIndex fineIndexAdj = m_fineLayout.findAdjacentIndex(a_finePatch, block, dir);
        addCoarseRegisterFromAdjacent(a_finePatch, a_coarseIndex, dir, fineIndexAdj);
    }
}

template<typename T, unsigned int C, MemType MEM>
void MBLevelFluxRegister<T,C,MEM>::addCoarseRegisterFromAdjacent(
    PatchID a_finePatch,
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
    int numFinePatches = m_fineLayout.localSize();
    m_fineRegisters.resize(numFinePatches);
    m_fineIndices.resize(numFinePatches);
    initFineRegisterIndices();
    for (auto fiter : m_fineLayout)
    {
        addFineRegisters(fiter);
    }
}

template<typename T, unsigned int C, MemType MEM>
void MBLevelFluxRegister<T,C,MEM>::initFineRegisterIndices ()
{
    for (auto& fineIndexArray : m_fineIndices)
    {
        for (int dd = 0; dd < DIM; dd++)
        {
            fineIndexArray[dd][(int)Side::Lo] = -1;
            fineIndexArray[dd][(int)Side::Hi] = -1;
        }
    }
}

template<typename T, unsigned int C, MemType MEM>
void MBLevelFluxRegister<T,C,MEM>::addFineRegisters(MBIndex a_fineIndex)
{
    PatchID finePatch = m_fineLayout.point(a_fineIndex);
    BlockIndex block = m_fineLayout.block(a_fineIndex);
    for (auto dir : Point::DirectionsOfCodim(1))
    {
        if (adjacentInProblemDomain(a_fineIndex, dir))
        {
            addFineRegisterFromAdjacent(a_fineIndex, dir);
        }
    }
}

template<typename T, unsigned int C, MemType MEM>
bool MBLevelFluxRegister<T,C,MEM>::adjacentInProblemDomain(
    MBIndex a_fineIndex,
    Point a_dir) const
{
    PatchID finePatch = m_fineLayout.point(a_fineIndex);
    BlockIndex block = m_fineLayout.block(a_fineIndex);
    Box patchDomain = m_fineLayout.patchDomain().getBlock(block).box();
    PatchID adjFinePatch = finePatch + a_dir;

    if (patchDomain.containsPoint(adjFinePatch)) { return true; }
    else if (m_fineLayout.onBlockBoundary(finePatch, block, a_dir)) { return true; }
    else { return false; }
}

template<typename T, unsigned int C, MemType MEM>
void MBLevelFluxRegister<T,C,MEM>::addFineRegisterFromAdjacent(
    MBIndex a_fineIndex,
    Point a_dir)
{
    int coord = a_dir.firstNonZeroIndex();
    int side = (a_dir[coord] < 0) ? 0 : 1;
    BlockIndex block = m_fineLayout.block(a_fineIndex);

    Point finePatch = m_fineLayout.point(a_fineIndex);
    MBIndex adjFineIndex = m_fineLayout.findAdjacentIndex(finePatch, block, a_dir);
    if (adjFineIndex == *m_fineLayout.end())
    {
        Point refRatio = m_refRatios[block];
        Box fineBox = m_fineLayout[a_fineIndex];
        Box coarseBox = fineBox.coarsen(refRatio);
        Box registerBox = coarseBox.adjacent(a_dir);

        shared_ptr<BoxData<T, C, MEM>> temp(new BoxData<T, C, MEM>(registerBox));
        Register<T, C, MEM> reg(temp, a_dir);

        int index = m_fineRegisters[a_fineIndex].size();
        m_fineRegisters[a_fineIndex].push_back(reg);
        m_fineIndices[a_fineIndex][coord][side] = index;
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
        int a_dir,
        const T& a_weight)
{
    PROTO_ASSERT(m_crseLayout.compatible(a_crseIndex),
        "MBLevelFluxRegister::incrementCoarseRegister | Error: Incompatible Data Index");
    BlockIndex block = m_crseLayout.block(a_crseIndex);
    auto& crseRegisters = m_crseRegisters[a_crseIndex];
    T weight = a_weight / m_dxCrse[block][a_dir];

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
        int a_dir,
        const T& a_weight)
{
    PROTO_ASSERT(m_fineLayout.compatible(a_fineIndex),
        "MBLevelFluxRegister::incrementCoarseRegister | Error: Incompatible Data Index");
    BlockIndex block = m_fineLayout.block(a_fineIndex);

    auto & fineRegister = m_fineRegisters[a_fineIndex];
    T weight = a_weight / m_dxCrse[block][a_dir];
    for (int k = 0; k < fineRegister.size(); k++)
    {
        auto& reg = fineRegister[k];
        if (reg.m_dir == a_dir)
        {
            auto& data = reg.m_data;
            auto  side = reg.m_side;
            (*data) += m_avgFlux[block][a_dir][side](a_flux, weight);
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
