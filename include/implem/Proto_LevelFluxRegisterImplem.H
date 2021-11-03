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
    std::cout << "Warning: FluxRegisterCopierOp::linearSize not implemented." << std::endl;
}
template<typename T, unsigned int C, MemType MEM>
void FluxRegisterCopierOp<T,C,MEM>::linearOut(
                                 void* a_buf,
                                 const Box& a_bx,
                                 const DataIndex& a_index) const
{
    // linearOut copies from source data holders into buffers.
    BoxData<T,C,MEM>& bd = m_register->sourceDB(a_bx,a_index);
    CInterval cint(0,C-1);
    bd.linearOut(a_buf, a_bx,cint);
}

template<typename T, unsigned int C, MemType MEM>
void FluxRegisterCopierOp<T,C,MEM>::linearIn(
                                 void* a_buf,
                                 const Box& a_bx,
                                 const DataIndex& a_index) const
{
    // linearIn copies from buffers into destination data holders.
    BoxData<T,C,MEM>& bd = m_register->destDB(a_bx,a_index);
    CInterval cint(0,C-1);
    bd.linearIn(a_buf, a_bx,cint);
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
    BoxData<T,C,MEM>& src = m_register->sourceBD(a_domain,a_domainIndex);
    BoxData<T,C,MEM>& dest =m_register->destBD(a_range,a_rangeIndex);
    Point shift = a_range.low() - a_domain.low();
    src.copyTo(dest, a_domain, shift);

#else
    // code for testing linearIn / linearOut while running serial.
    PR_TIMERS("LOP_localCopy_2");
    int sizeSource = size(a_domain);
    int sizeDest   = size(a_range);

    if (sizeSource != sizeDest)
    {
        MayDay<void>::Abort("LinearizationTest failure: dest and source have different sizes");
    }

    vector<char> buffer(sizeSource);
    {
        PR_TIMERS("linearIn / linearOut inside LDOP::op");
        void* b = (void*)&(buffer[0]);
        linearOut(b, a_domain,a_domainIndex);
        linearIn(b, a_range,a_rangeIndex);
    }
#endif
}

template<typename T, unsigned int C, MemType MEM>
void FluxRegisterCopier<T,C,MEM>::buildMotionPlans (FluxRegisterCopierOp<T, C, MEM >& a_op)
{
    auto & lfr = *(a_op.m_register);
    auto & refRatio = lfr.m_refRatio;
    auto & dblFine = lfr.m_dblFine;
    auto & dblCoarse = lfr.m_dblCoarse;
    auto & coarseRegisters = lfr.m_coarseRegisters;
    auto & hasFineRegister = lfr.m_hasFineRegister;
    auto & fineRegisters = lfr.m_fineRegisters;

    DataIterator ddit(dblCoarse);
    for (ddit.begin();*ddit != ddit.end();++ddit)
    {
        int myindex = dblCoarse.myIndex(*ddit);
        auto & cReg = coarseRegisters[myindex];
        for (int k = 0; k < cReg.size(); k++)
        {
            auto & reg = cReg[k];
            // Find fine index corresponding to this flux register.
            Box destBox = reg.m_data->box();
            int lohi = reg.m_lohi;
            int dir = reg.m_dir;
            Point ptFine = destBox.low()*refRatio/dblFine.boxSize() + Point::Basis(dir,1-2*lohi);
            DataIndex diFine = dblFine.find(ptFine);
            PR_assert(diFine != dblFine.end());
            // Compute fine register box (nontrivial b/c of possibility of periodic BCs. 
            Box coarsenedFineBox = dblFine.box(diFine).coarsen(refRatio);
            Box srcBox;
            if (lohi == 0)
            {
                srcBox = coarsenedFineBox.adjCellLo(dir,1);
            }
            else
            {
                srcBox =  coarsenedFineBox.adjCellHi(dir,1);
            }                       
            MotionItem* toMotionItem = new MotionItem(diFine, *ddit, srcBox, destBox);
            auto fineID = dblFine.procID(diFine);
            if (fineID == procID())
            { 
                this->m_localMotionPlan.push_back(toMotionItem);
            }
            else
            {
                toMotionItem->procID =fineID;
                this->m_toMotionPlan.push_back(toMotionItem);
            }
        }
    }
    // Iterate over fine grids that I own to obtain entries in m_fromMotionPlan.           
    DataIterator sdit(dblFine);
    int destProcID = procID();
    for (sdit.begin();*sdit != sdit.end();++sdit)
    {
        int myindex = dblFine.myIntIndex(*sdit);
        Box coarsenedFineBox = dblFine.box(*sdit).coarsen(refRatio);
        for (int dir = 0; dir < DIM; dir++)
        {
            for (int lohi = 0;lohi < 2; lohi++)
            {
                if (hasFineRegister[myindex][dir][lohi] >= 0)
                {
                    Box srcBox;
                    // compute box for possible register.
                    if (lohi == 0)
                    {
                        srcBox = coarsenedFineBox.adjCellLo(dir,1);
                    }
                    else
                    {
                        srcBox = coarsenedFineBox.adjCellHi(dir,1);
                    }
                    auto ddit = dblCoarse.find(srcBox.low()/dblCoarse.boxSize());
                    int coarseProcID = dblCoarse.procID(ddit);
                    Point shiftCoarse =
                        dblCoarse.problemDomain().shifted(srcBox.low()) - srcBox.low();
                    Box destBox = srcBox.shift(shiftCoarse);
                    if (coarseProcID != procID())
                    {
                        MotionItem* fromMotionItem = 
                            new MotionItem(*sdit,ddit,srcBox,destBox);
                        fromMotionItem->procID = coarseProcID;
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
        const DisjointBoxLayout& a_dblFine,
        const DisjointBoxLayout& a_dblCoarse,
        const Point& a_refRatio)
{
    define(a_dblFine,a_dblCoarse,a_refRatio);
}

template<typename T, unsigned int C, MemType MEM>
void LevelFluxRegister<T,C,MEM>::define (
        const DisjointBoxLayout& a_dblFine,
        const DisjointBoxLayout& a_dblCoarse,
        const Point& a_refRatio)
{
    m_refRatio = a_refRatio;
    m_dblFine = a_dblFine;
    m_dblCoarse = a_dblCoarse;


    // Set up Coarse data holders.
    DataIterator ditC(m_dblCoarse);

    Point finePatchesBoxSize = m_dblCoarse.boxSize()*m_refRatio / m_dblFine.boxSize();
    // Check to see whether Coarse patch is exactly covered by fine patches.
    PR_assert(
            finePatchesBoxSize*m_dblFine.boxSize() == m_dblCoarse.boxSize()*m_refRatio);
    // Check to see whether coarsened fine patch is more than one coarse cell wide in all directions that are being refined. 

    // Fine registers are defined relative to the fine patch coordinates, and coarse 
    // registers are defined relative to the coarse patch coordinates. What happens in 
    // the presence of periodic boundaries is mediated in LevelFluxRegister<T,C,MEM>Operator.

    for (int dir = 0; dir < DIM; dir++)
    {
        if (m_refRatio[dir] > 1)
        {
            PR_assert(m_dblFine.boxSize()[dir]/m_refRatio[dir] > 1);
        }
    }
    m_coarseRegisters.clear();
    m_destData.clear();
    m_coarseIndices.clear();
    m_fineRegisters.clear();
    m_hasFineRegister.clear();

    int mySizeC = m_dblCoarse.mySize();
    m_coarseRegisters.resize(mySizeC);
    m_destData.resize(mySizeC);
    m_coarseIndices.resize(mySizeC);
    // Set up coarse data holders.  
    for (ditC.begin();*ditC!= ditC.end();++ditC)
    {
        int myindex = m_dblCoarse.myIndex(*ditC);
        Box bxC = m_dblCoarse[*ditC];
        int k = 0;
        m_coarseRegisters[myindex].clear();
        m_destData[myindex].clear();
        m_coarseIndices[myindex].clear();
        //setup: compute box corresponding to possible fine patches covered by coarse patch.
        Box finePatchPoints = bxC.refine(m_refRatio).coarsen(m_dblFine.boxSize());
        BoxIterator bit(finePatchPoints);
        // Iterate over the fine patches contained in the coarse box, and see whether patches 
        // adjacent to them in each direction can contribute to a flux register in the coarse 
        // patch. 
        for (bit.begin();bit.done();++bit)
        {
            Point pt = *bit;
            DataIndex di = m_dblFine.index(pt);
            Box cfBox,destBox;
            for (int d = 0; d < DIM; d++)
            {
                for (int lohi = 0; lohi < 2;lohi++)
                {
                    Point ptAdj = *bit + Point::Basis(d,2*lohi-1);
                    DataIndex diAdj = m_dblFine.find(ptAdj);
                    if (diAdj != m_dblFine.end() && di == m_dblFine.end())
                    {
                        Box cfBox = 
                            m_dblFine[di].coarsen(m_refRatio);
                        if (lohi == 0)
                        {
                            destBox = cfBox.adjCellHi(d,1);
                        }
                        else
                        {
                            destBox = cfBox.adjCellLo(d,1);
                        }
                    }
                    // Intersect with the Problem domain to shift the Box to its proper 
                    // periodic image in the coarse patch.

                    destBox = m_dblCoarse.problemDomain()&destBox; 

                    // Create register and insert it into the coarse register data structure.

                    int thisKey = key(cfBox,di);
                    shared_ptr<BoxData<T,C> > temp1(new BoxData<T,C>(destBox));
                    shared_ptr<BoxData<T,C> > temp2(new BoxData<T,C>(destBox));
                    Register<T,C,MEM> reg(temp1,d,1);
                    m_coarseRegisters[myindex].push_back(reg);
                    m_destData[myindex].push_back(temp2);
                    m_coarseIndices[myindex][thisKey] = k;
                    k++;
                }
            }
        }
    }  

    // Set up fine data holders.
    DataIterator ditF(m_dblFine);
    int mySizeF = m_dblFine.mySize();
    m_fineRegisters.resize(mySizeF);
    m_hasFineRegister.resize(mySizeF);
    for (ditF.begin();*ditF!=ditF.end();++ditF)
    {
        Point pt = m_dblFine.point(*ditF);
        Box bx = Box(pt,pt).refine(m_dblFine.boxSize());
        int myIndex = m_dblFine.myIndex(*ditF);
        int k = 0;
        Box srcBox;
        for (int d = 0; d < DIM; d++)
        {
            for (int lohi= 0;lohi < 2; lohi++)
            {
                Point srcPt;
                if (lohi == 0) 
                {
                    srcPt = pt + Point::Basis(d,-1);
                    srcBox = bx.adjCellLo(d,1);
                }
                else 
                {
                    srcPt = pt + Point::Basis(d);
                    srcBox = bx.adjCellHi(d,1);
                }
                DataIndex adjIndex = m_dblFine.index(srcPt);
                if (adjIndex == m_dblFine.end())
                {
                    shared_ptr<BoxData<T,C> > temp(new BoxData<T,C>(srcBox)); 
                    Register<T,C,MEM> reg(temp,d,lohi);
                    m_fineRegisters[myIndex].push_back(reg);
                    m_hasFineRegister[myIndex][d][lohi] = k;
                    k++;
                }
            }
        }
    }

    // We now have enoughinformation to build the motion plans for the copier.
    FluxRegisterCopierOp<T,C,MEM> op(*this);
    m_copier.define(op);
    m_copier.buildMotionPlans(op);

    // Create stencils for averaging down fluxes.
    for (int d = 0; d < DIM; d++)
    {
        for (int lohi = 0; lohi < 1; lohi++)
        {
            Stencil<double> & sten = m_avgFlux[d][lohi] = (0.0)*Shift(Point::Zeros());
            Box bx(Point::Zeros(),(Point::Ones() - Point::Basis(d))*m_refRatio);
            BoxIterator bit(bx);
            T coef = 1.0/bx.size();
            for (auto iter = bx.begin(); iter !=bx.end(); ++iter)

            {
                sten += coef*Shift(*iter);
            }
            sten.srcRatio() = m_refRatio; 
            if (lohi == 0)
            {
                sten.destShift() = Point::Basis(d,-1);
            }
        }
    }
}

template<typename T, unsigned int C, MemType MEM>
void LevelFluxRegister<T,C,MEM>::incrementCoarse (
        const BoxData<T,C,MEM>& a_flux,
        const DataIndex& a_di,
        const T& a_weight,
        unsigned int a_dir)
{
    std::cout << "FluxRegister::incrementCoarse: dir = " << a_dir << std::endl;
    std::cout << "m_coarseRegisters.size(): " << m_coarseRegisters.size() << std::endl; 
    int myindex = m_dblCoarse.myIndex(a_di);
    auto & cReg = m_coarseRegisters[myindex];
    std::cout << "m_coarseRegisters[" << myindex << "].size(): " << cReg.size() << std::endl;
    //Stencil<double> scaleLow  = (-a_weight)*Shift(Point::Zeros() - Point::Basis(a_dir));
    //Stencil<double> scaleHigh =   a_weight*Shift(Point::Zeros());
    Stencil<double> scaleLow = -a_weight*Shift::Basis(a_dir, -1);
    Stencil<double> scaleHigh = a_weight*Shift::Zeros();
    for (int k = 0; k < cReg.size();k++)
    {   
        if (cReg[k].m_dir == a_dir)
        {
            auto & reg = cReg[k].m_data; 
            int side = cReg[k].m_lohi;
            if (side == 0)
            {
                (*reg) += scaleLow(a_flux);
            }
            else
            {
                (*reg) += scaleHigh(a_flux);
            }
        }
    }
}
template<typename T, unsigned int C, MemType MEM>
void LevelFluxRegister<T,C,MEM>::incrementFine (
        const BoxData<T, C , MEM>& a_flux,
        const DataIndex& a_di,
        const T& a_weight,
        unsigned int a_dir)
{
    int myIndex = m_dblFine.myIndex(a_di);
    auto & myRegs = m_fineRegisters[myIndex];
    for (int k = 0; k < myRegs.size(); k++)
    {
        const auto & rg = myRegs[k];
        auto dir = rg.m_dir;
        if (dir == a_dir)
        {
            auto lohi = rg.m_lohi;
            auto db = rg.m_data;
            (*db) += m_avgFlux[dir][lohi](a_flux);
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

    // Data has been copied to m_destData. We compute the difference between the 
    // coarse and fine registers, and increment a_coarseData.

    Stencil<double> increment = (-a_weight)*Shift(Point::Zeros());
    for (auto dit = a_coarseData.begin();*dit != dit.end();++dit)
    {
        int ind = m_dblCoarse.myIntIndex(*dit);
        int size = m_destData[ind].size();
        auto & dest = m_destData[ind];
        auto & creg = m_coarseRegisters[ind];
        for (int k = 0; k < size; k++)
        {
            *(dest[k]) += increment(*(creg[k].m_data));  
        }
        for (int k = 0; k < size; k++)
        {
            a_coarseData[*dit] += *(creg[k].m_data);
        }
    }          
}

template<typename T, unsigned int C, MemType MEM>
void LevelFluxRegister<T,C,MEM>::reset()
{
    // set values in the registers to zero.
    for (auto dit = m_dblCoarse.begin();*dit != dit.end(); ++dit)
    {
        //int myIndex = dit.myIndex();
        int myIndex = *dit;
        for (int k = 0; k < m_coarseIndices[myIndex].size();k++)
        {
            m_coarseRegisters[myIndex][k].m_data->setToZero();
            m_destData[myIndex][k]->setToZero();
        }
    }
    for (auto dit = m_dblFine.begin();*dit != dit.end(); ++dit)
    {
        //int myIndex = dit.myIndex();
        int myIndex = *dit;
        for (int k = 0; k < m_fineRegisters[myIndex].size();k++)
        {
            m_fineRegisters[myIndex][k].m_data->setToZero();
        }
    }
}

    template<typename T, unsigned int C, MemType MEM>
int LevelFluxRegister<T,C,MEM>::key(const Box& a_bx,const DataIndex& a_di)
{
    // Find which direction is normal to the register. We use the fact that the register
    // is of size 1 only in the normal direction.

    int normal = -1;
    for (int dir = 0; dir < DIM; dir++)
    {
        if (a_bx.size(dir) == 1) normal = dir;
    }
    PR_assert(normal !=-1);
    Point coarsePt = m_dblCoarse.point(a_di);
    Box bxCoarse(Point::Zeros(),m_dblCoarse.boxSize());
    Point shiftedPt = a_bx.low() - coarsePt;
    unsigned int ind = bxCoarse.index(shiftedPt);
    unsigned int sz = bxCoarse.size();

    // Find whether this is a low or high side register. We use the fact that, if 
    // fine grids are at least two coarse cells wide, there can be only one fine patch 
    // per direction adjacent to the register Box.

    Point ptFinePatch = a_bx.low()*m_refRatio/m_dblFine.boxSize();
    Box finePatchPoints = Box(ptFinePatch,ptFinePatch).
        refine(m_dblCoarse.boxSize()).refine(m_refRatio).coarsen(m_dblFine.boxSize());
    Point ptLow = ptFinePatch - Point::Basis(normal);
    Point ptHigh = ptFinePatch + Point::Basis(normal);
    DataIndex diLow = m_dblFine.find(ptLow);
    PR_assert(diLow != m_dblFine.end() ||  m_dblFine.find(ptHigh) != m_dblFine.end());
    int lohi = 1;
    if (diLow != m_dblFine.end()) lohi = 0;

    return ind + sz*normal + sz*DIM*lohi;
}

template<typename T, unsigned int C, MemType MEM>
BoxData<T,C,MEM>& LevelFluxRegister<T,C,MEM>::sourceBD(const Box& a_bx,const DataIndex& a_di)
{
    Point ptFine = m_dblFine.point(a_di);
    Box bxfine = Box(ptFine,ptFine).refine(m_dblFine.boxSize());
    int myindex = m_dblFine.myIntIndex(a_di);
    int lohi = -1;
    int dir = -1;
    for (int d = 0; d< DIM;d++)
    {
        if (a_bx == bxfine.adjCellLo(d,1))
        {
            lohi = 0;
            dir = d;
        }
        else if (a_bx == bxfine.adjCellHi(d,1))
        {
            lohi = 1;
            dir = d;
        }
    }
    PR_assert(dir != -1);
    int k = m_hasFineRegister[dir][lohi][myindex];
    return *(m_fineRegisters[myindex][k].m_data);
}

template<typename T, unsigned int C, MemType MEM>
BoxData<T,C,MEM>& LevelFluxRegister<T,C,MEM>::destBD(const Box& a_bx,const DataIndex& a_di)
{
    int keyCoarse = key(a_bx,a_di);
    int myindex = m_dblCoarse.myIndex(a_di);
    return *(m_coarseRegisters[myindex][keyCoarse].m_data);
}

template<typename T, unsigned int C, MemType MEM>
void LevelFluxRegister<T,C,MEM>::print(std::string a_name) const
{
    std::cout << "Level Flux Register: " << a_name << endl;
}
