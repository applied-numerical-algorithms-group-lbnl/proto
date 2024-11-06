template<typename T, unsigned int C, MemType MEM>
MBFluxRegisterCopierOp<T,C,MEM>::MBFluxRegisterCopierOp(MBLevelFluxRegister<T, C, MEM>& a_register)
{
    m_register = &a_register;
}

template<typename T, unsigned int C, MemType MEM>
int MBFluxRegisterCopierOp<T,C,MEM>::linearSize(
                                   const Box& a_bx,
                                   const DataIndex<MBBoxPartition>& a_index) const
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

    PR_TIMERS("LOP_localCopy_1");
    auto& src = m_register->sourceData(domain, srcIndex);
    auto& dst = m_register->destData(  range,  dstIndex);
    Point shift = range.low() - domain.low();
    src.copyTo(dst, domain, shift);
#else
    // code for testing linearIn / linearOut while running serial.
    PR_TIMERS("LOP_localCopy_2");
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
    auto & refRatio = lfr.m_refRatio;
    auto & fineLayout = lfr.m_fineLayout;
    auto & crseLayout = lfr.m_crseLayout;
    auto & crseRegisters = lfr.m_crseRegisters;
    auto & crseIncrement = lfr.m_crseIncrement;
    auto & fineRegisters = lfr.m_fineRegisters;
    auto & fineIndices = lfr.m_fineIndices;

    //std::cout << "Building Destination Motion Items" << std::endl;
    for (auto citer = crseLayout.begin(); citer.ok(); ++citer)
    {
        auto& localRegisters = crseRegisters[*citer];
        for (int k = 0; k < localRegisters.size(); k++)
        {
            // Find the index of the fineRegister source.
            // The index is associated with the fine patch contributing the averaged flux
            auto& reg     = localRegisters[k];
            Box   destBox = reg.m_data->box();
            auto  side    = reg.m_side;
            int   dir     = reg.m_dir;
            Point finePatchPoint = destBox.low()*refRatio/fineLayout.boxSize() + Point::Basis(dir,side);
            //std::cout << "\tCoarse Box: " << destBox << " | Dir: " << dir << " | Side: " << side << std::endl;
            //std::cout << "\tAdjacent Fine Patch: " << finePatchPoint << std::endl;
            auto fineIndex = fineLayout.find(finePatchPoint);
            
            PROTO_ASSERT(fineIndex != *fineLayout.end(),
                "MBFluxRegisterCopier::buildMotionPlans | Error: Could not find fine register index.");
            // Compute fine register box (nontrivial b/c of possibility of periodic BCs. 
            Box cfBox  = fineLayout.box(fineIndex).coarsen(refRatio);
            Box srcBox = cfBox.adjacent(dir, flip(side), 1);
            // Boxes for the fine registers are defined relative to the fine patch
            // coordinates, even if the register Boxes are periodic images outside the
            // ProblemDomain Box - PC.
            MBMotionItem toMotionItem(fineIndex, *citer, srcBox, destBox);
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
    // std::cout << "Building Source Motion Items" << std::endl;
    for (auto fiter = fineLayout.begin(); fiter.ok(); ++fiter)
    {
      //Box cfBox = fiter.box().coarsen(refRatio);
        for (int dir = 0; dir < DIM; dir++)
        {
            SideIterator siter;
            for (siter.begin(); siter.ok(); ++siter)
            {
                int k = fineIndices[*fiter][dir][*siter];
                if (k >= 0) //otherwise there is no fine register for this (dir, side)
                {
                  Box srcBox = fineRegisters[*fiter][k].m_data->box();
                  Point crsePatchPoint = srcBox.low()/crseLayout.boxSize();
                    // std::cout << "\tCoarse Box: " << srcBox << " | Dir: " << dir << " | Side: " << *siter << std::endl;
                    
                    auto crseIndex  = crseLayout.find(crsePatchPoint);
                    int  crseProcID = crseLayout.procID(crseIndex);
                    //Automatically defines dstBox to be the correct periodic image. NO IT DOESN'T.
                    Box dstBox = crseLayout.domain() & srcBox;
                    dstBox = crseLayout.domain().image(dstBox); 
                    if (crseProcID != procID())
                    {
                        MBMotionItem fromMotionItem(*fiter, crseIndex, srcBox, dstBox);
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
        const MBDisjointBoxLayout& a_fineLayout,
        const MBDisjointBoxLayout& a_crseLayout,
        const std::vector<Point>&   a_refRatios,
        std::vector<Array<T, DIM>>& a_dxCoarse)
{
    define(a_fineLayout,a_crseLayout,a_refRatio,a_dxCrse);
}

template<typename T, unsigned int C, MemType MEM>
MBLevelFluxRegister<T,C,MEM>::MBLevelFluxRegister(
        const MBDisjointBoxLayout& a_fineLayout,
        const MBDisjointBoxLayout& a_crseLayout,
        Point&          a_refRatio,
        Array<T, DIM>&  a_dxCoarse)
{
    define(a_fineLayout, a_crseLayout, a_refRatio, a_dxCoarse);
}

template<typename T, unsigned int C, MemType MEM>
void MBLevelFluxRegister<T,C,MEM>::define (
        const MBDisjointBoxLayout& a_crseLayout,
        const MBDisjointBoxLayout& a_fineLayout,
        Point&          a_refRatio,
        Array<T, DIM>&  a_dxCoarse)
{
    int numBlocks = a_crseLayout.numBlocks();
    std::vector<Point> refRatios(a_refRatio, numBlocks);
    std::vector<Array<T, DIM>> dx(a_dxCoarse, numBlocks);
    define(a_fineLayout, a_crseLayout, refRatios, dx);
}

template<typename T, unsigned int C, MemType MEM>
void MBLevelFluxRegister<T,C,MEM>::define (
        const MBDisjointBoxLayout& a_crseLayout,
        const MBDisjointBoxLayout& a_fineLayout,
        const std::vector<Point>&   a_refRatios,
        std::vector<Array<T, DIM>>& a_dxCoarse)
{
    int numBlocks = a_crseLayout.numBlocks();
    PROTO_ASSERT(a_fineLayout.numBlocks() == numBlocks,
            "MBLevelFluxRegister::define | Error: layouts have different block counts.");
    PROTO_ASSERT(a_refRatios.size() == numBlocks,
            "MBLevelFluxRegister::define | Error: wrong number of refinement ratios.");
    PROTO_ASSERT(a_dxCoarse.size() == numBlocks,
            "MBLevelFluxRegister::define | Error: wrong number of grid spacings.");

    // Fine registers are defined relative to the fine patch coordinates, and coarse 
    // registers are defined relative to the coarse patch coordinates. What happens in 
    // the presence of periodic boundaries is mediated in MBLevelFluxRegister<T,C,MEM>Operator.
    
    m_refRatio = a_refRatio;
    m_dxCrse = a_dxCrse;
    m_fineLayout = a_fineLayout;
    m_crseLayout = a_crseLayout;

    // Check to see whether Coarse patch is exactly covered by fine patches.
    for (int bi = 0; bi < numBlocks; bi++)
    {
        Point ri = m_refRatio[bi];
        Point patchRefRatio = 
            m_crseLayout.boxSize(bi)*ri / m_fineLayout.boxSize(bi);
        PROTO_ASSERT(patchRefRatio*m_fineLayout.boxSize(bi) == m_crseLayout.boxSize(bi)*ri,
            "MBLevelFluxRegister::define | Error: coarse region is not covered by fine region.");
    }

    // Check to see whether coarsened fine patch is more than one coarse cell wide in all directions that are being refined. 
    for (int bi = 0; bi < numBlocks; bi++)
    {
        for (int dir = 0; dir < DIM; dir++)
        {
            Point ri = m_refRatio[bi];
            if (ri[dir] > 1)
            {
                PROTO_ASSERT(m_fineLayout.boxSize(bi)[dir]/ri[dir] > 1,
                    "MBLevelFluxRegister::define | Error:\
                    fine patch cell size is smaller than the refinement ratio.");
            }
        }
    }
    m_crseRegisters.clear();
    m_crseIncrement.clear();
    m_crseIndices.clear();
    m_fineRegisters.clear();
    m_fineIndices.clear();

    int numCrsePatches = m_crseLayout.localSize();
    m_crseRegisters.resize(numCrsePatches);
    m_crseIncrement.resize(numCrsePatches);
    m_crseIndices.resize(numCrsePatches);
    
    auto& crsePatchDomain = m_crseLayout.partition().patchDomain();
    auto& finePatchDomain = m_fineLayout.partition().patchDomain();

    for (auto citer : m_crseLayout)
    {
        BlockIndex block = m_crseLayout.block(citer);
        Point refRatio = m_refRatio[block];
        Point fineBoxSize = m_fineLayout.boxSize(block);
        Point crseBoxSize = m_crseLayout.boxSize(block);
        Box crseBox = m_crseLayout[citer];
        m_crseRegisters[citer].clear();
        m_crseIncrement[citer].clear();
        m_crseIndices[citer].clear();

        Box finePatchPoints = crseBox.refine(refRatio).coarsen(fineBoxSize);
        PROTO_ASSERT(finePatchDomain.contains(finePatchPoints),
                "MBLevelFluxRegister::define | Error: Data corruption.");
        int k = 0;
        for (auto finePoint : finePatchPoints)
        {
            // check if fine patch is ABSENT from fine layout
            auto fineIndex = m_fineLayout.find(pt, block);
            if (fineIndex != *m_fineLayout.end()) { continue; }
        
            // check adjacent fine patch that is PRESENT in fine layout
            Box cfBox = Box(finePoint, finePoint).refine(fineBoxSize).coarsen(refRatio);
            for (auto dir : Box::Kernel(1))
            {
                if (dir.codim() != 1) { continue; } //flux register ignores codim > 1 boundaries
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
                    continue;
                } else {
                    // INVALID case
                    PROTO_ASSERT(false,
                        "MBLevelFluxRegister::define | Error: Data corruption.");
                }
                if (fineIndexAdj == *m_fineLayout.end()) { continue; } //adj patch is absent

                // found a valid cf boundary in block interior
                Box destBox = cfBox.face(dir);

                int thisKey = key(destBox, citer);
                shared_ptr<BoxData<T,C,MEM> > temp1(new BoxData<T,C,MEM>(destBox));
                shared_ptr<BoxData<T,C,MEM> > temp2(new BoxData<T,C,MEM>(destBox));
                Register<T,C,MEM> reg(temp1, dir);
                m_crseRegisters[citer].push_back(reg);
                m_crseIncrement[citer].push_back(temp2);
                m_crseIndices  [citer][thisKey] = k;
                k++;
            } // end dir loop
        } // end fine patch loop
    } // end coarse loop

    // Set up fine data holders.
    int numFinePatches = m_fineLayout.localSize();
    m_fineRegisters.resize(numFinePatches);
    m_fineIndices.resize(numFinePatches);
    for (auto fiter : m_fineLayout)
    {
        BlockIndex block = m_fineLayout.block(fiter);
        Point finePoint = m_fineLayout.point(fiter);
        Box fineBox = m_fineLayout[fiter];
        Box cfBox = fineBox.coarsen(m_refRatio);
        int k = 0;
        for (auto dd = 0; dd < DIM; dd++)
        {
            SideIterator siter;
            for (siter.begin(); siter.ok(); ++siter())
            {
                Point dir = Point::Basis(dd,*siter);
                Point finePointAdj = finePoint + dir;
                auto fineIndexAdj = *m_fineLayout.end();
                if (finePatchDomain.inInterior(finePointAdj, block))
                {
                    // INTERIOR case (normal case)
                    fineIndexAdj = m_fineLayout.find(finePointAdj, block);
                } else if (finePatchDomain.inBlockBoundary(srcPoint, block))
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
                if (adjIndex == *m_fineLayout.end())
                {
                    //found a cf boundary between finePoint and finePointAdj
                    Box srcBox = cfBox.adjacent(dir);

                    shared_ptr<BoxData<T,C,MEM> > temp(new BoxData<T,C,MEM>(srcBox));
                    Register<T,C,MEM> reg(temp, d, *siter);
                    m_fineRegisters[fiter].push_back(reg);
                    m_fineIndices[fiter][d][(int)(*siter)] = k;
                    k++;
                } else {
                    m_fineIndices[*fiter][d][(int)(*siter)] = -1;
                }

            }
        }
    }
#if 0
    for (auto citer : m_crseLayout)
    {
        BlockIndex block = m_crseLayout.block(citer);
        Point refRatio = m_refRatio[block];
        Point fineBoxSize = m_fineLayout.boxSize(block);
        Point crseBoxSize = m_crseLayout.boxSize(block);
        Box crseBox = m_crseLayout[citer];
        m_crseRegisters[citer].clear();
        m_crseIncrement[citer].clear();
        m_crseIndices[citer].clear();
        

        //setup: compute box corresponding to possible fine patches covered by coarse patch.
        Box finePatchPoints = crseBox.refine(refRatio).coarsen(fineBoxSize);
        // Iterate over the fine patches contained in the coarse box, and see whether patches 
        // adjacent to them in each direction can contribute to a flux register in the coarse 
        // patch. 
        
        int k = 0;
        for (auto pt : finePatchPoints)
        {
            auto fineIndex = m_fineLayout.find(pt);
            // only check fine patches NOT in part of the fine layout
            if (fineIndex != *m_fineLayout.end()) { continue; } 
            Box cfBox = Box(pt, pt).refine(fineBoxSize).coarsen(refRatio);
            for (int dir = 0; dir < DIM; dir++)
            {
                SideIterator siter;
                for (siter.begin(); siter.ok(); ++siter)
                {
                    Point ptAdj = pt + Point::Basis(dir,*siter);
                    auto adjIndex = m_fineLayout.find(ptAdj);
                    if (adjIndex != *m_fineLayout.end())
                    {
                        // fineIndex is the index of a fine patch NOT in the fine layout
                        // adjIndex is the index of a fine patch IN the fine layout
                        // hence, the interface between the two is a coarse-fine boundary
                        Box destBox = cfBox.face(dir, *siter);
                        // Intersect with the Problem domain to shift the Box to its proper 
                        // periodic image in the coarse patch.
                        destBox = m_crseLayout.domain() & destBox; 
#if PR_VERBOSE > 1
                        pr_out() << "\t\t\t\tCreating Register: destBox = " << destBox << std::endl;
#endif
                        //std::cout << "\t\t\tFound Boundary Patch: " << ptAdj << " | Shifted Coarsened Box: " << destBox << std::endl;

                        // Create register and insert it into the coarse register data structure.
                        // the register is not uniquely determined by it's (dir, side) 
                        // the key is generated using destBox which incodes (dir, side, tile position)
                        int thisKey = key(destBox, *citer);
                        shared_ptr<BoxData<T,C,MEM> > temp1(new BoxData<T,C,MEM>(destBox));
                        shared_ptr<BoxData<T,C,MEM> > temp2(new BoxData<T,C,MEM>(destBox));
                        Register<T,C,MEM> reg(temp1, dir, *siter);
                        m_crseRegisters[*citer].push_back(reg);
                        m_crseIncrement[*citer].push_back(temp2);
                        m_crseIndices  [*citer][thisKey] = k;
                        k++;
                    }
                } // end for side
            } // end for DIM
        } // end for fine in coarse
    }  // end for coarse
    // Set up fine data holders.
    int numFinePatches = m_fineLayout.localSize();
    m_fineRegisters.resize(numFinePatches);
    m_fineIndices.resize(numFinePatches);
    for (auto fiter = m_fineLayout.begin(); fiter.ok(); ++fiter)
    {
        Point finePatchPoint = m_fineLayout.point(*fiter);
        Box fineBox = m_fineLayout[*fiter];
        Box cfBox = fineBox.coarsen(m_refRatio);
        int k = 0;
        for (int d = 0; d < DIM; d++)
        {
            SideIterator siter;
            for (siter.begin(); siter.ok(); ++siter)
            {
                Point srcPt = finePatchPoint + Point::Basis(d, *siter);
                auto adjIndex = m_fineLayout.find(srcPt);
                bool srcPatchMissing = (adjIndex == *m_fineLayout.end());
                bool srcPatchInDomain = m_fineLayout.domain().contains(srcPt);
                // look for adjacent patches NOT in the fine layout

                if (srcPatchMissing && srcPatchInDomain)
                {
                    // Fine registers are AVERAGED fine fluxes
                    // hence, domains are coarse
                    Box srcBox = cfBox.adjacent(d, *siter, 1);
                    // Deal with possible periodic images
                    srcBox = m_crseLayout.domain() & srcBox;
                    shared_ptr<BoxData<T,C,MEM> > temp(new BoxData<T,C,MEM>(srcBox));
                    Register<T,C,MEM> reg(temp, d, *siter);
                    m_fineRegisters[*fiter].push_back(reg);
                    m_fineIndices[*fiter][d][(int)(*siter)] = k;
                    k++;
                } else {
                    m_fineIndices[*fiter][d][(int)(*siter)] = -1;
                }
            } // end for side
        } // end for DIM
    } // end for fine

#endif
    // We now have enough information to build the motion plans for the copier.
    MBFluxRegisterCopierOp<T,C,MEM> op(*this);
    m_copier.define(op);

    // Create stencils for averaging down fluxes.
    for (int d = 0; d < DIM; d++)
    {   
        SideIterator siter;
        for (siter.begin(); siter.ok(); ++siter)
        {
            // Side::Lo is not a typo -CLG
            auto avg = Stencil<T>::AvgDownFace(d, Side::Lo, m_refRatio);
            if (*siter == Side::Hi)
            {
                avg = avg*(-1);
            } else if (*siter == Side::Lo)
            {
                // Side::Lo is not a typo here either -CLG
                avg.destShift() = Point::Basis(d, Side::Lo);
            }
            m_avgFlux[d][*siter] = avg;
        }
    }

    reset();
    #if 0
    for (auto iter = m_crseLayout.begin(); iter.ok(); ++iter)
    {
        int numRegisters = m_crseRegisters[*iter].size();
        for (int ii = 0; ii < numRegisters; ii++)
        {
            auto& reg = m_crseRegisters[*iter][ii];
        }
    }
    for (auto iter = m_fineLayout.begin(); iter.ok(); ++iter)
    {
        int numRegisters = m_fineRegisters[*iter].size();
        for (int ii = 0; ii < m_fineRegisters[*iter].size(); ii++)
        {
            auto& reg = m_fineRegisters[*iter][ii];
        }
    }
#endif
}

template<typename T, unsigned int C, MemType MEM>
void MBLevelFluxRegister<T,C,MEM>::incrementCoarse (
        const BoxData<T,C,MEM>& a_flux,
        const DataIndex<MBBoxPartition>& a_crseIndex,
        const T& a_weight,
        unsigned int a_dir)
{
    //std::cout << "FluxRegister::incrementCoarse | patch index: " << (int)a_crseIndex << " | box: " << a_flux.box() << std::endl;
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
void MBLevelFluxRegister<T,C,MEM>::incrementFine (
        const BoxData<T, C , MEM>& a_flux,
        const DataIndex<MBBoxPartition>& a_fineIndex,
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
void MBLevelFluxRegister<T,C,MEM>::reflux(
        LevelBoxData<T,C,MEM>& a_coarseData,
        const T& a_weight)
{
    PROTO_ASSERT(m_crseLayout.compatible(a_coarseData.layout()),
            "MBLevelFluxRegister::reflux | Error: \
            Input data has incompatible layout.");
    
    m_copier.execute();

    // fine register data has been copied to m_crseIncrement.
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

template<typename T, unsigned int C, MemType MEM>
void MBLevelFluxRegister<T,C,MEM>::reset(T a_val)
{
    // set values in the registers to zero.
    for (auto citer = m_crseLayout.begin(); citer.ok(); ++citer)
    {
        for (int k = 0; k < m_crseIndices[*citer].size(); k++)
        {
            m_crseRegisters[*citer][k].m_data->setVal(a_val);
            m_crseIncrement[*citer][k]->setToZero();
        }
    }
    for (auto fiter = m_fineLayout.begin(); fiter.ok(); ++fiter)
    {
        for (int k = 0; k < m_fineRegisters[*fiter].size();k++)
        {
            m_fineRegisters[*fiter][k].m_data->setToZero();
        }
    }
}

// a_bx is the coarse box updated by refluxing
template<typename T, unsigned int C, MemType MEM>
int MBLevelFluxRegister<T,C,MEM>::key(const Box& a_bx, const DataIndex<MBBoxPartition>& a_di)
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

    /* 
    int dir = -1;
    int side = -1;
    Box crseBox = m_crseLayout[a_di];
    for (int d = 0; d < DIM; d++)
    {
        SideIterator siter;
        for (siter.begin(); siter.ok(); ++siter)
        {
            Box faceBox = crseBox.face(d, *siter, 1);
            if (faceBox.contains(a_bx))
            {
                Point tileSize = a_bx.high() - a_bx.low() + Point::Ones() - Point::Basis(d);
                Point shiftPoint = faceBox.low();
                Box tile = a_bx.shift(-shiftPoint).coarsen(tileSize);
                Box tileDomain = faceBox.shift(-shiftPoint).coarsen(tileSize);
                PROTO_ASSERT(tileDomain.contains(tile), "MBLevelFluxRegister::key | Error: Could not compute tile index");
                PROTO_ASSERT(tileDomain.contains(tile), "MBLevelFluxRegister::key | Error: Invalid tile box");
                int numTiles = tileDomain.size();
                int tileIndex = tileDomain.index(tile.low());

                return tileIndex + numTiles*d + numTiles*DIM*(*siter);
            }
        }
    }
    MayDay<void>::Abort("MBLevelFluxRegister::key | Error: could not generate key.");
    */
}

template<typename T, unsigned int C, MemType MEM>
BoxData<T,C,MEM>& MBLevelFluxRegister<T,C,MEM>::sourceData(const Box& a_bx, const DataIndex<MBBoxPartition>& a_fineIndex)
{
    std::shared_ptr<BoxData<T,C,MEM>> retval;
    Box fineBox = m_fineLayout.box(a_fineIndex);
    Box cfBox = fineBox.coarsen(m_refRatio);
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
BoxData<T,C,MEM>& MBLevelFluxRegister<T,C,MEM>::destData(const Box& a_bx, const DataIndex<MBBoxPartition>& a_crseIndex)
{
    int crseKey       = key(a_bx, a_crseIndex);
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
