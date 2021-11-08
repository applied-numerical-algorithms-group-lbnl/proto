template<typename T, unsigned int C, MemType MEM>
FluxRegisterCopierOp<T,C,MEM>::FluxRegisterCopierOp(LevelFluxRegister<T, C, MEM>& a_register)
{
    m_register = &a_register;
}

template<typename T, unsigned int C, MemType MEM>
int FluxRegisterCopierOp<T,C,MEM>::linearSize(
                                   const Box& a_bx,
                                   const DataIndex& a_index) const
{
    return C*a_bx.size()*sizeof(T);
}
template<typename T, unsigned int C, MemType MEM>
void FluxRegisterCopierOp<T,C,MEM>::linearOut(
                                 void* a_buf,
                                 const Box& a_bx,
                                 const DataIndex& a_index) const
{
    // linearOut copies from source data holders into buffers.
    BoxData<T,C,MEM>& data = m_register->sourceData(a_bx, a_index);
    CInterval cint(0,C-1);
    data.linearOut(a_buf, a_bx, cint);
}

template<typename T, unsigned int C, MemType MEM>
void FluxRegisterCopierOp<T,C,MEM>::linearIn(
                                 void* a_buf,
                                 const Box& a_bx,
                                 const DataIndex& a_index) const
{
    // linearIn copies from buffers into destination data holders.
    BoxData<T,C,MEM>& data = m_register->destData(a_bx,a_index);
    CInterval cint(0,C-1);
    data.linearIn(a_buf, a_bx,cint);
}

template<typename T, unsigned int C, MemType MEM>
void FluxRegisterCopierOp<T,C,MEM>::localCopy(
        const Box& a_domain,
        const DataIndex& a_domainIndex,
        const Box& a_range,
        const DataIndex& a_rangeIndex) const
{
#define noLinearizationTest 1

#if noLinearizationTest

    PR_TIMERS("LOP_localCopy_1");
    BoxData<T,C,MEM>& src  = m_register->sourceData(a_domain,a_domainIndex);
    BoxData<T,C,MEM>& dest = m_register->destData(a_range,a_rangeIndex);
    Point shift = a_range.low() - a_domain.low();
    src.copyTo(dest, a_domain, shift);

#else
    // code for testing linearIn / linearOut while running serial.
    PR_TIMERS("LOP_localCopy_2");
    int srcSize = linearSize(a_domain, a_domainIndex);
    int dstSize = linearSize(a_range,  a_rangeIndex);

    if (srcSize != dstSize)
    {
        MayDay<void>::Abort("LinearizationTest failure: dest and source have different sizes");
    }

    vector<char> buffer(sizeSource);
    {
        PR_TIMERS("linearIn / linearOut inside LDOP::op");
        void* b = (void*)&(buffer[0]);
        linearOut(b, a_domain, a_domainIndex);
        linearIn (b, a_range,  a_rangeIndex);
    }
#endif
}

template<typename T, unsigned int C, MemType MEM>
void FluxRegisterCopier<T,C,MEM>::buildMotionPlans (FluxRegisterCopierOp<T, C, MEM >& a_op)
{
    auto & lfr = *(a_op.m_register);
    auto & refRatio = lfr.m_refRatio;
    auto & fineLayout = lfr.m_fineLayout;
    auto & crseLayout = lfr.m_crseLayout;
    auto & crseRegisters = lfr.m_crseRegisters;
    auto & crseIncrement = lfr.m_crseIncrement;
    auto & fineRegisters = lfr.m_fineRegisters;
    auto & fineIndices = lfr.m_fineIndices;

    for (auto citer = crseLayout.begin(); citer.ok(); ++citer)
    {
        auto& localRegisters = crseRegisters[*citer];
        auto& localIncrement = *(crseIncrement[*citer]);
        for (int k = 0; k < localRegisters.size(); k++)
        {
            // Find the index of the fineRegister source.
            // The index is associated with the fine patch contributing the averaged flux
            auto& reg     = localRegisters[k];
            Box   destBox = reg.m_data->box();
            auto  side    = reg.m_side;
            int   dir     = reg.m_dir;
            Point finePatchPoint = destBox.low()*refRatio/fineLayout.boxSize() + Point::Basis(dir,side);
            DataIndex fineIndex = fineLayout.find(finePatchPoint);
            PROTO_ASSERT(fineIndex != fineLayout.end(),
                "FluxRegisterCopier::buildMotionPlans | Error: Could not find fine register index.");
            
            // Compute fine register box (nontrivial b/c of possibility of periodic BCs. 
            //Box cfBox  = fineLayout.box(fineIndex).coarsen(refRatio);
            //Box srcBox = cfBox.adjacent(dir, flip(side), 1);
            Box srcBox = destBox; // define uses exclusively shifted boxes so periodicity is already taken care of -CLG
            MotionItem* toMotionItem = new MotionItem(fineIndex, *citer, srcBox, destBox);
            auto fineID = fineLayout.procID(fineIndex);
            if (fineID == procID())
            { 
                this->m_localMotionPlan.push_back(toMotionItem);
            }
            else
            {
                toMotionItem->procID = fineID;
                this->m_toMotionPlan.push_back(toMotionItem);
            }
        } // for each local coarse register
    } // end for coarse patch

    // Iterate over fine grids that I own to obtain entries in m_fromMotionPlan.           
    int destProcID = procID();
    for (auto fiter = fineLayout.begin(); fiter.ok(); ++fiter)
    {
        Box cfBox = fiter.box().coarsen(refRatio);
        for (int dir = 0; dir < DIM; dir++)
        {
            SideIterator siter;
            for (siter.begin(); siter.ok(); ++siter)
            {
                int k = fineIndices[*fiter][dir][*siter];
                if (k >= 0) //otherwise there is no fine register for this (dir, side)
                {
                    Box srcBox = cfBox.adjacent(dir, *siter, 1);
                    // have to account for periodicity since we called adjacent
                    srcBox = crseLayout.domain() & srcBox;
                    Point crsePatchPoint = srcBox.low()/crseLayout.boxSize();
                    auto crseIndex  = crseLayout.find(crsePatchPoint);
                    int  crseProcID = crseLayout.procID(crseIndex);
                    Box dstBox = srcBox; // just for clarity
                    if (crseProcID != procID())
                    {
                        MotionItem* fromMotionItem = new MotionItem(*fiter, crseIndex, srcBox, destBox);
                        fromMotionItem->procID = crseProcID;
                        this->m_fromMotionPlan.push_back(fromMotionItem);
                    }
                }
            }
        }
    }     
    this->sort();
}

template<typename T, unsigned int C, MemType MEM>
LevelFluxRegister<T,C,MEM>::LevelFluxRegister(
        const DisjointBoxLayout& a_fineLayout,
        const DisjointBoxLayout& a_crseLayout,
        const Point& a_refRatio)
{
    define(a_fineLayout,a_crseLayout,a_refRatio);
}

template<typename T, unsigned int C, MemType MEM>
void LevelFluxRegister<T,C,MEM>::define (
        const DisjointBoxLayout& a_fineLayout,
        const DisjointBoxLayout& a_crseLayout,
        const Point& a_refRatio)
{
    // Fine registers are defined relative to the fine patch coordinates, and coarse 
    // registers are defined relative to the coarse patch coordinates. What happens in 
    // the presence of periodic boundaries is mediated in LevelFluxRegister<T,C,MEM>Operator.
    
    m_refRatio = a_refRatio;
    m_fineLayout = a_fineLayout;
    m_crseLayout = a_crseLayout;

    // Check to see whether Coarse patch is exactly covered by fine patches.
    Point patchRefRatio = m_crseLayout.boxSize()*m_refRatio / m_fineLayout.boxSize();
    PROTO_ASSERT(patchRefRatio*m_fineLayout.boxSize() == m_crseLayout.boxSize()*m_refRatio
        "LevelFluxRegister::define | Error: coarse region is not covered exactly by fine region.");

    // Check to see whether coarsened fine patch is more than one coarse cell wide in all directions that are being refined. 
    for (int dir = 0; dir < DIM; dir++)
    {
        if (m_refRatio[dir] > 1)
        {
            PROTO_ASSERT(m_fineLayout.boxSize()[dir]/m_refRatio[dir] > 1
                "LevelFluxRegister::define | Error: fine patch cell size is smaller than the refinement ratio.");
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
    
    // Set up coarse data holders.  
    for (auto citer = m_crseLayout.begin(); citer.ok(); ++citer)
    {
        Box crseBox = citer.box();
        m_crseRegisters[*citer].clear();
        m_crseIncrement[*citer].clear();
        m_crseIndices[*citer].clear();
        
        //setup: compute box corresponding to possible fine patches covered by coarse patch.
        Box finePatchPoints = crseBox.refine(m_refRatio).coarsen(m_fineLayout.boxSize());
        // Iterate over the fine patches contained in the coarse box, and see whether patches 
        // adjacent to them in each direction can contribute to a flux register in the coarse 
        // patch. 
        int k = 0;
        for (auto biter = finePatchPoints.begin(); biter.ok(); ++biter)
        {
            Point pt = *biter;
            DataIndex fineIndex = m_fineLayout.index(pt);
            // only check fine patches NOT in part of the fine layout
            if (fineIndex != m_fineLayout.end()) { continue; } 
            Box cfBox = Box(pt, pt).refine(m_fineLayout.boxSize).coarsen(m_refRatio);
            for (int d = 0; d < DIM; d++)
            {
                SideIterator siter;
                for (siter.begin(); siter.ok(); ++siter)
                {
                    Point ptAdj = *biter + Point::Basis(d,*siter);
                    // find accounts for periodic boundaries
                    DataIndex adjIndex = m_fineLayout.find(ptAdj);
                    if (adjIndex != m_fineLayout.end())
                    {
                        // fineIndex is the index of a fine patch NOT in the fine layout
                        // adjIndex is the index of a fine patch IN the fine layout
                        // hence, the interface between the two is a coarse-fine boundary
                        Box destBox = cfBox.edge(dir, *siter);
                        // Intersect with the Problem domain to shift the Box to its proper 
                        // periodic image in the coarse patch.
                        destBox = m_crseLayout.problemDomain() & destBox; 

                        // Create register and insert it into the coarse register data structure.
                        // the register is not uniquely determined by it's (dir, side) 
                        // the key is generated using destBox which incodes (dir, side, tile position)
                        int thisKey = key(destBox, *citer);
                        shared_ptr<BoxData<T,C> > temp1(new BoxData<T,C>(destBox));
                        shared_ptr<BoxData<T,C> > temp2(new BoxData<T,C>(destBox));
                        Register<T,C,MEM> reg(temp1, d, *siter);
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
        Point finePatchPoint = fiter.point();
        Box fineBox = fiter.box();
        int k = 0;
        for (int d = 0; d < DIM; d++)
        {
            SideIterator siter;
            for (siter.begin(); siter.ok(); ++siter)
            {
                Point srcPt = finePatchPoint + Point::Basis(d, *siter);
                DataIndex adjIndex = m_fineLayout.index(srcPt);
                // look for adjacent patches NOT in the fine layout
                if (adjIndex == m_fineLayout.end())
                {
                    // Fine registers are AVERAGED fine fluxes
                    // hence, domains are coarse
                    Box srcBox = fineBox.coarsen(m_refRatio);
                    srcBox = srcBox.adjacent(d, *siter, 1);
                    // Deal with possible periodic images
                    srcBox = m_crseLayout.domain() & srcBox;
                    shared_ptr<BoxData<T,C> > temp(new BoxData<T,C>(srcBox));
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

    // We now have enough information to build the motion plans for the copier.
    FluxRegisterCopierOp<T,C,MEM> op(*this);
    m_copier.define(op);
    // m_copier.buildMotionPlans(op);

    // Create stencils for averaging down fluxes.
    for (int d = 0; d < DIM; d++)
    {   
        SideIterator siter;
        for (siter.begin(); siter.ok(); ++siter)
        {
            auto avg = Stencil<T>::AvgDownEdge(d, m_refRatio);
            if (*siter == Side::Lo)
            {
                //FIXME: Is this correct?
                avg.destShift() = Point::Basis(d, -1);
            }
            m_avgFlux[d][(int)(*siter)] = avg;
        }
    }
}

template<typename T, unsigned int C, MemType MEM>
void LevelFluxRegister<T,C,MEM>::incrementCoarse (
        const BoxData<T,C,MEM>& a_flux,
        const DataIndex& a_crseIndex,
        const T& a_weight,
        unsigned int a_dir)
{
    auto& crseRegister = m_crseRegisters[a_crseIndex];
    Stencil<double> scaleLow = -a_weight*Shift::Basis(a_dir, -1);
    Stencil<double> scaleHigh = a_weight*Shift::Zeros();
    for (int k = 0; k < crseRegister.size(); k++)
    {   
        auto& reg = crseRegisters[k];
        if (reg.m_dir == a_dir)
        {
            auto& data = reg.m_data; 
            auto  side = reg.m_side;
            // FIXME: Make sure registers are set to zero when defined
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
void LevelFluxRegister<T,C,MEM>::incrementFine (
        const BoxData<T, C , MEM>& a_flux,
        const DataIndex& a_fineIndex,
        const T& a_weight,
        unsigned int a_dir)
{
    auto & fineRegister = m_fineRegisters[a_fineIndex];
    for (int k = 0; k < fineRegister.size(); k++)
    {
        auto& reg = fineRegister[k];
        if (reg.m_dir == a_dir)
        {
            auto& data = reg.m_data;
            auto  side = reg.m_side
            (*data) += m_avgFlux[a_dir][(int)(side)](a_flux);
        }
    }
}

template<typename T, unsigned int C, MemType MEM>
void LevelFluxRegister<T,C,MEM>::reflux(
        LevelBoxData<T,C,MEM>& a_coarseData,
        const T& a_weight)
{
    // invoke Copier.
    m_copier.execute();

    // Data has been copied to m_crseIncrement. We compute the difference between the 
    // coarse and fine registers, and increment a_coarseData.
    Stencil<double> increment = (-a_weight)*Shift::Zeros();
    for (auto iter = a_coarseData.begin(); iter.ok(); ++iter)
    {
        int numRegisters = m_crseIncrement[*iter].size();
        auto& dest = m_crseIncrement[*iter];
        auto& creg = m_crseRegisters[*iter];
        // subtract off the coarse fluxes
        for (int k = 0; k < numRegisters; k++)
        {
            *(dest[k]) += increment(*(creg[k].m_data));  
        }
        // add increment to input data
        for (int k = 0; k < numRegisters; k++)
        {
            a_coarseData[*iter] += *(dest[k]);
        }
    }          
}

template<typename T, unsigned int C, MemType MEM>
void LevelFluxRegister<T,C,MEM>::reset()
{
    // set values in the registers to zero.
    for (auto citer = m_crseLayout.begin(); citer.ok(); ++citer)
    {
        for (int k = 0; k < m_crseIndices[*citer].size(); k++)
        {
            m_crseRegisters[*citer][k].m_data->setToZero();
            m_crseIncrement[*citer][k]->setToZero();
        }
    }
    for (auto fiter = m_fineLayout.begin(); fiter.ok(); ++fiter)
    {
        for (int k = 0; k < m_fineRegisters[myIndex].size();k++)
        {
            m_fineRegisters[*fiter][k].m_data->setToZero();
        }
    }
}

// a_bx is the coarse box updated by refluxing
template<typename T, unsigned int C, MemType MEM>
int LevelFluxRegister<T,C,MEM>::key(const Box& a_bx, const DataIndex& a_di)
{
    // Find which direction is normal to the register. We use the fact that the register
    // is of size 1 only in the normal direction.

    int normal = -1;
    for (int dir = 0; dir < DIM; dir++)
    {
        if (a_bx.size(dir) == 1) { normal = dir; }
    }
    PROTO_ASSERT(normal != -1,
        "LevelFluxRegister::key | Error: Could not determine normal direction of input Box.");
    Point coarsePt = m_crseLayout.point(a_di);
    Box bxCoarse(Point::Zeros(),m_crseLayout.boxSize());
    Point shiftedPt = a_bx.low() - coarsePt;
    unsigned int ind = bxCoarse.index(shiftedPt);
    unsigned int sz = bxCoarse.size();

    // Find whether this is a low or high side register. We use the fact that, if 
    // fine grids are at least two coarse cells wide, there can be only one fine patch 
    // per direction adjacent to the register Box.

    Point ptFinePatch = a_bx.low()*m_refRatio/m_fineLayout.boxSize();
    Box finePatchPoints = Box(ptFinePatch,ptFinePatch).
        refine(m_crseLayout.boxSize()).refine(m_refRatio).coarsen(m_fineLayout.boxSize());
    Point ptLow = ptFinePatch - Point::Basis(normal);
    Point ptHigh = ptFinePatch + Point::Basis(normal);
    DataIndex diLow = m_fineLayout.find(ptLow);
    PR_assert(diLow != m_fineLayout.end() ||  m_fineLayout.find(ptHigh) != m_fineLayout.end());
    int lohi = 1;
    if (diLow != m_fineLayout.end()) lohi = 0;

    return ind + sz*normal + sz*DIM*lohi;
}

template<typename T, unsigned int C, MemType MEM>
BoxData<T,C,MEM>& LevelFluxRegister<T,C,MEM>::sourceData(const Box& a_bx, const DataIndex& a_fineIndex)
{
    Box fineBox = m_fineLayout.box(a_fineIndex);
    int fineIntIndex = m_fineLayout.localIntIndex(a_fineIndex);
    SideIterator siter;
    for (siter.begin(); siter.ok(); ++siter)
    {
        for (int dir = 0; dir < DIM; dir++)
        {
            if (a_bx == fineBox.adjacent(dir, *siter, 1))
            {
                int k = m_fineIndices[dir][(int)(*siter)][fineIntIndex];
                return *(m_fineRegisters[fineIntIndex][k].m_data);
            }
        }
    }
    std::cout << "Invalid Box: " << a_bx << " or invalid index: " << fineIntIndex << std::endl;
    PROTO_ASSERT(false,
        "LevelFluxRegister::sourceData | Error: Could not locate data.");
}

template<typename T, unsigned int C, MemType MEM>
BoxData<T,C,MEM>& LevelFluxRegister<T,C,MEM>::destData(const Box& a_bx, const DataIndex& a_crseIndex)
{
    int crseKey       = key(a_bx, a_crseIndex);
    int localIntIndex = m_crseLayout.localIntIndex(a_crseIndex);
    int mapIndex      = m_crseIndices[localIntIndex][crseKey];
    return *(m_crseIncrement[localIntIndex][mapIndex]);
}

template<typename T, unsigned int C, MemType MEM>
void LevelFluxRegister<T,C,MEM>::print(std::string a_name) const
{
    std::cout << "Level Flux Register: " << a_name << endl;
}
