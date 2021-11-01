template<typename T, unsigned int C, MemType MEM, unsigned int CTR>
LevelBoxData<T, C, MEM, CTR>::
LevelBoxData()
{
    m_isDefined = false;
}

template<typename T, unsigned int C, MemType MEM, unsigned int CTR>
LevelBoxData<T, C, MEM, CTR>::
LevelBoxData(const DisjointBoxLayout& a_layout, 
             const Point& a_ghost)
{
    define(a_layout, a_ghost);
}
template<typename T, unsigned int C, MemType MEM, unsigned int CTR>
void LevelBoxData<T, C, MEM, CTR>::define(
                 const DisjointBoxLayout& a_layout, 
                 const Point& a_ghost)
{
    m_isDefined = true;
    m_ghost = a_ghost;
    m_dbl = a_layout;
    for (auto dit = begin(); dit.ok(); ++dit)
    {
        //We don't allow LevelBoxData to use StackAllocator in allocating BoxDatas
        bool use_stackalloc = false;
        auto patch = std::make_shared<BoxData<T,C,MEM>>(dit.box().grow(m_ghost), use_stackalloc);
        m_data.push_back(patch);
    }

    if (m_ghost != Point::Zeros())
    {
        m_exchangeCopier.define(LevelCopierOp<T, C, MEM, CTR>(*this, *this));
    }
}

template<typename T, unsigned int C, MemType MEM, unsigned int CTR>
template<typename Func, typename... Srcs>
void LevelBoxData<T, C, MEM, CTR>::initialize(Func& a_func, Srcs... a_srcs)
{
    for (auto iter = this->begin(); iter.ok(); ++iter)
    {
        auto& patch = this->operator[](*iter);
        forallInPlace_p(a_func, patch, a_srcs...);
    }
}

template<typename T, unsigned int C, MemType MEM, unsigned int CTR>
unsigned int
LevelBoxData<T, C, MEM, CTR>::patchSize() const
{
    int size = 1;
    Point boxSize = layout().boxSize() + 2*ghost();
    for (int ii = 0; ii < DIM; ii++)
    {
        size *= boxSize[ii];
    }
    return C*size;
}

template<typename T, unsigned int C, MemType MEM, unsigned int CTR>
unsigned int
LevelBoxData<T, C, MEM, CTR>::offset(int a_proc) const
{
    Point boxSizeVect = layout().boxSize() + 2*ghost();
    int boxSize = 1;
    for (int ii = 0; ii < DIM; ii++)
    {
        boxSize *= boxSizeVect[ii];
    }
    return layout().offset(a_proc)*boxSize*C;
}

template<typename T, unsigned int C, MemType MEM, unsigned int CTR>
unsigned int 
LevelBoxData<T, C, MEM, CTR>::size() const
{
    return m_dbl.size();
}

template<typename T, unsigned int C, MemType MEM, unsigned int CTR>
void   
LevelBoxData<T, C, MEM, CTR>::setToZero()
{
    for (auto iter = begin(); iter.ok(); ++iter)
    {
        auto& patch = this->operator[](*iter);
        patch.setVal(0);
    }
}

template<typename T, unsigned int C, MemType MEM, unsigned int CTR>
void 
LevelBoxData<T, C, MEM, CTR>::exchange()
{
    PR_TIME("LevelBoxData::exchange");
#ifdef PR_MPI
    {
        PR_TIME("MPI_Barrier exchange");
        barrier();
    }
#endif
    m_exchangeCopier.execute();
}

template<typename T, unsigned int C, MemType MEM, unsigned int CTR>
void 
LevelBoxData<T, C, MEM, CTR>::
copyTo(LevelBoxData<T, C, MEM, CTR>& a_dest) const
{
    // we trust that copier will keep *this constant by only moving data to dst
    auto mutableSrc = const_cast<LevelBoxData<T, C, MEM, CTR>*> (this);
    LevelCopierOp<T, C, MEM, CTR> op(*mutableSrc, a_dest);
    LevelCopier<T, C, MEM, CTR> copier;
    copier.define(op);
    copier.execute();
}

template<typename T, unsigned int C, MemType MEM, unsigned int CTR>
void 
LevelBoxData<T, C, MEM, CTR>::linearIn(
    void* a_buffer,
    const Box& a_box,
    const DataIndex& a_index)
{
    auto& patch = this->operator[](a_index);
    CInterval comps(0, C-1);
    patch.linearIn(a_buffer, a_box, comps);
}

template<typename T, unsigned int C, MemType MEM, unsigned int CTR>
void 
LevelBoxData<T, C, MEM, CTR>::linearIn(
    void* a_buffer,
    const DataIndex& a_index)
{
    Box b = this->operator[](a_index).box();
    linearIn(a_buffer, b, a_index);
}

template<typename T, unsigned int C, MemType MEM, unsigned int CTR>
void 
LevelBoxData<T, C, MEM, CTR>::linearIn( void* a_buffer)
{
    char* readPoint = (char*)a_buffer;
    for (auto iter = this->begin(); *iter != iter.end(); ++iter)
    {
        linearIn(readPoint, *iter);
        readPoint += linearSize(*iter); 
    }
}

template<typename T, unsigned int C, MemType MEM, unsigned int CTR>
void 
LevelBoxData<T, C, MEM, CTR>::linearOut(
    void* a_buffer,
    const Box& a_box,
    const DataIndex& a_index) const
{
    auto& patch = this->operator[](a_index);
    CInterval comps(0, C-1);
    patch.linearOut(a_buffer, a_box, comps);
}

template<typename T, unsigned int C, MemType MEM, unsigned int CTR>
void 
LevelBoxData<T, C, MEM, CTR>::linearOut(
    void* a_buffer,
    const DataIndex& a_index) const
{
    Box b = this->operator[](a_index).box();
    linearOut(a_buffer, b, a_index);
}

template<typename T, unsigned int C, MemType MEM, unsigned int CTR>
void 
LevelBoxData<T, C, MEM, CTR>::linearOut( void* a_buffer) const
{
    char* writePoint = (char*)a_buffer;
    for (auto iter = this->begin(); *iter != iter.end(); ++iter)
    {
        linearOut(writePoint, *iter);
        writePoint += linearSize(*iter); 
    }
}

template<typename T, unsigned int C, MemType MEM, unsigned int CTR>
size_t
LevelBoxData<T, C, MEM, CTR>::linearSize(
    const Box& a_box,
    const DataIndex& a_index) const
{
    return C*a_box.size()*sizeof(T);
}

template<typename T, unsigned int C, MemType MEM, unsigned int CTR>
size_t
LevelBoxData<T, C, MEM, CTR>::linearSize( const DataIndex& a_index) const
{
    Box b = this->operator[](a_index).box();
    return this->linearSize(b, a_index);
}

template<typename T, unsigned int C, MemType MEM, unsigned int CTR>
size_t
LevelBoxData<T, C, MEM, CTR>::linearSize() const
{
    size_t totalSize = 0;
    for (auto iter = this->begin(); *iter != iter.end(); ++iter)
    {
        totalSize += size(*iter);
    }
    return totalSize;
}

//// Interaction with iterators.
template<typename T, unsigned int C, MemType MEM, unsigned int CTR>
DataIterator LevelBoxData<T, C, MEM, CTR>::begin() const
{
    DataIterator retval(m_dbl);
    retval.begin();
    return retval;
} 

template<typename T, unsigned int C, MemType MEM, unsigned int CTR>
BoxData<T,C,MEM>& 
LevelBoxData<T, C, MEM, CTR>::operator[](const DataIndex& a_di) 
{
    return (*(m_data[m_dbl.myIntIndex(a_di)]));
}

template<typename T, unsigned int C, MemType MEM, unsigned int CTR>
const BoxData<T,C,MEM>& 
LevelBoxData<T, C, MEM, CTR>::operator[](const DataIndex& a_di) const
{
    return (*(m_data[m_dbl.myIntIndex(a_di)]));
}

// =======================================================================
// LEVEL COPIER OP

template<typename T, unsigned int C, MemType MEM, unsigned int CTR>
LevelCopierOp<T, C, MEM, CTR>::LevelCopierOp(
        LevelBoxData<T, C, MEM, CTR>& a_src,
        LevelBoxData<T, C, MEM, CTR>& a_dst)
{
    m_src = &a_src;
    m_dst = &a_dst;
}

template<typename T, unsigned int C, MemType MEM, unsigned int CTR>
int
LevelCopierOp<T, C, MEM, CTR>::linearSize(const Box& a_box,
        const DataIndex& a_index) const
{
    //TODO: Does this do the right thing?
    return C*a_box.size()*sizeof(T);
}

template<typename T, unsigned int C, MemType MEM, unsigned int CTR>
void
LevelCopierOp<T, C, MEM, CTR>::linearOut(void* a_buf, 
        const Box& a_bx,
        const DataIndex& a_index) const
{
    auto& patch = m_src->operator[](a_index);
    CInterval cint(0,C-1);
    patch.linearOut(a_buf, a_bx, cint);
}

template<typename T, unsigned int C, MemType MEM, unsigned int CTR>
void
LevelCopierOp<T, C, MEM, CTR>::linearIn(void* a_buf, 
        const Box& a_bx,
        const DataIndex& a_index) const
{
    auto& patch = m_dst->operator[](a_index);
    CInterval cint(0,C-1);
    patch.linearIn(a_buf, a_bx, cint);
}

template<typename T, unsigned int C, MemType MEM, unsigned int CTR>
void
LevelCopierOp<T, C, MEM, CTR>::localCopy(
        const Box& a_domain, const DataIndex& a_domainIndex,
        const Box& a_range,  const DataIndex& a_rangeIndex) const
{
    int LinearizationTest = 1;
    if (LinearizationTest == 0)
    {
        PR_TIMERS("LOP_localCopy_1");
        auto& src = m_src->operator[](a_domainIndex);
        auto& dst = m_dst->operator[](a_rangeIndex);
        Point shift = a_range.low() - a_domain.low();
        src.copyTo(dst, a_domain, shift);
    }
    else
    {
        PR_TIMERS("LOP_localCopy_2");
        int sizeSrc = linearSize(a_domain, a_domainIndex);
        int sizeDst = linearSize(a_range, a_rangeIndex);

        if (sizeSrc != sizeDst)
        {
            MayDay<void>::Abort("LinearizationTest failure: dest and source have different sizes");
        }

        vector<char> buffer(sizeSrc);
        {
            PR_TIMERS("linearIn / linearOut inside LDOP::op");
            void* b = (void*)&(buffer[0]);
            linearOut(b, a_domain, a_domainIndex);
            linearIn (b, a_range,  a_rangeIndex);
        }
    }
}

// =======================================================================
// LEVEL COPIERS
template<typename T, unsigned int C, MemType MEM, unsigned int CTR>
void
LevelCopier<T, C, MEM, CTR>::buildMotionPlans(LevelCopierOp<T, C, MEM, CTR>& a_op)
{
    // Loop over boxes in the destination to construct toMotionPlan, localMotionPlan
    int destProcID = procID();                                                              // get dest procID
    for (auto diter = a_op.m_dst->begin(); diter.ok(); ++diter)                             // for each patch in dest
    {
        Point ghost = a_op.m_dst->ghost();
        auto bx = diter.box().grow(ghost);                                                  // get dest layout box
        NeighborIterator nit((*a_op.m_src).layout(), bx);                                   // make neighbor iterator from src layout and dest box
        for (nit.begin(); nit.ok(); ++nit)                                                  // for each neighbor of dest box
        {
            int nbrProcID = nit.procID();                                                   // get neighbor procID
            Box destBx = nit.destBox();                                                     // 
            if (!destBx.empty())
            {
                MotionItem* toMotionItem = new MotionItem(*nit, *diter, destBx, destBx);
                if (nbrProcID == procID())
                {
                    // local move
                    this->m_localMotionPlan.push_back(toMotionItem);
                } else {
                    toMotionItem->procID = nbrProcID;
                    this->m_toMotionPlan.push_back(toMotionItem);
                }
            }
        } // end neighbor loop
    } // end dst loop
    
    for (auto siter = a_op.m_src->begin(); siter.ok(); ++siter)
    {
        Box bx = siter.box();
        NeighborIterator nit((*a_op.m_dst).layout(), bx);
        for (nit.begin(); nit.ok(); ++nit)
        {
            int nbrProcID = nit.procID();
            Box srcBx = nit.srcBox();
            if (nbrProcID != destProcID)
            {
                MotionItem* fromMotionItem = new MotionItem(*siter, *nit, srcBx, srcBx);
                fromMotionItem->procID = nbrProcID;
                this->m_fromMotionPlan.push_back(fromMotionItem);
            }
        }
    } // end src loop
    this->sort();  
}

template<typename T, unsigned int C, MemType MEM, unsigned int CTR>
void
LevelExchangeCopier<T, C, MEM, CTR>::buildMotionPlans(LevelCopierOp<T, C, MEM, CTR>& a_op)
{
    auto layout = a_op.m_src->layout();
    auto ghost  = a_op.m_src->ghost();

    int myProcID = Proto::procID();
    for (auto dit = layout.begin(); dit.ok(); ++dit)
    {
        Box bx0 = dit.box();
        Box bx = layout.problemDomain() & dit.box().grow(ghost);
        NeighborIterator nit(layout, bx);
        for (nit.begin(); nit.ok(); ++nit)
        {
            int nbrProcID = nit.procID();
            Point nbrShift = nit.shift();
            // For each pair of interacting boxes, we construct a
            // Motion to, from, and local.
            // Compute the boxes corresponding to filling from *dit to *nit
            // The source box is shifted, the dest box is not.
            
            Box destBx = nit.destBox();
            if (!destBx.empty() && !bx0.contains(destBx))
            {
                Box srcBx = nit.srcBox();
                MotionItem* toMotionItem = new MotionItem(*nit, *dit, srcBx, destBx);
                if (nbrProcID == myProcID)
                {
                    // local motion
                    this->m_localMotionPlan.push_back(toMotionItem);
                } else {
                    toMotionItem->procID = nbrProcID;
                    this->m_toMotionPlan.push_back(toMotionItem);
                }
            }
            
            // Construct the motion for filling the ghost region for *nit from *dit

            if (nbrProcID != myProcID)
            {
                Box nbrbx = Box(nit.point(), nit.point()).refine(layout.boxSize()).grow(ghost);
                Box nbrbxShift = nbrbx.shift(-nbrShift);
                Box srcFromBx = nbrbxShift & bx0;
                Box destFromBx = srcFromBx.shift(nbrShift);
                MotionItem* fromMotionItem
                    = new MotionItem(*dit, *nit, srcFromBx, destFromBx);
                fromMotionItem->procID = nbrProcID;
                this->m_fromMotionPlan.push_back(fromMotionItem);
            }
        } // end neighbor loop
    } // end patch loop
    this->sort();
}

template<typename Func, typename... Srcs>
void forallInPlace(Func& a_F, Srcs&&... a_srcs)
{

}



/* FIXME: This is all still here in case I need to reference the previous (working) implementation
template<typename T, unsigned int C, MemType MEM, unsigned int CTR>
void 
LevelBoxData<T, C, MEM, CTR>::makeItSo(
              const Copier&       a_copier,
              const LDOP&         a_op) const
{
  if(s_verbosity > 0) 
    {
      pout() << "makeit so copier = " << endl;
      //a_copier.print();
    }
  {
    PR_TIME("makeItSoBegin");
    makeItSoBegin(a_copier, a_op);
  }
  {
    PR_TIME("makeItSoLocalCopy");
    makeItSoLocalCopy(a_copier, a_op);
  }
  {
    PR_TIME("makeItSoEnd");
    makeItSoEnd(a_op);
  }
}
template<typename T, unsigned int C, MemType MEM, unsigned int CTR>
void LevelBoxData<T, C, MEM, CTR>::makeItSoBegin(
                        const Copier&       a_copier,
                        const LDOP&         a_op) const
{
  // The following five functions are nullOps in uniprocessor mode

#ifdef PR_MPI
  allocateBuffers(a_copier,a_op);  //monkey with buffers, set up 'fromMe' and 'toMe' queues

  writeSendDataFromMeIntoBuffers(a_op);

  // If there is nothing to recv/send, don't go into these functions
  // and allocate memory that will not be freed later.  (ndk)
  // The #ifdef PR_MPI is for the m_buff->m_toMe and m_buff->m_fromMe
  {
    PR_TIME("post_messages");
    this->m_buff->numReceives = m_buff->m_toMe.size();

    if (this->m_buff->numReceives > 0)
      {
        postReceivesToMe(); // all non-blocking
      }  

    this->m_buff->numSends = m_buff->m_fromMe.size();
    if (this->m_buff->numSends > 0)
      {
        postSendsFromMe();  // all non-blocking
      }
  }    
#endif 
}
template<typename T, unsigned int C, MemType MEM, unsigned int CTR>
void 
LevelBoxData<T, C, MEM, CTR>::makeItSoLocalCopy(const Copier&       a_copier,
                       const LDOP&         a_op) const
{
  PR_TIME("local copying");
  CopyIterator it(a_copier, CopyIterator::LOCAL);  
  int items=it.size();
  for (int n=0; n<items; n++)
    {
      const MotionItem& item = it[n];
      //debugging bit to force serial code to run parallel bits
#if 0  
      size_t bufsize_src = a_op.size(item.fromRegion);
      size_t bufsize_dst = a_op.size(item.fromRegion);
      if(bufsize_src != bufsize_dst)
        {
          MayDay<void*>::Error("buffer size mismatch");
        }
      char* charbuffer = new char[bufsize_src];
      a_op.linearOut(charbuffer, item.fromRegion,item.fromIndex);
      a_op.linearIn (charbuffer, item.toRegion,item.toIndex);
        
      delete[] charbuffer;
#else 
      a_op.op(item.fromRegion,item.fromIndex,item.toRegion,item.toIndex);
#endif
        
    }
}
template<typename T, unsigned int C, MemType MEM, unsigned int CTR>
void 
LevelBoxData<T, C, MEM, CTR>::makeItSoEnd(const LDOP& a_op) const
{
#ifdef PR_MPI
  // Uncomment and Move this out of unpackReceivesToMe()  (ndk)
  //cout << "completePendingSends" << " , proc = " << procID() << endl;
  completePendingSends(); // wait for sends from possible previous operation
  //cout << "unpackReceivesToMe"  << " , proc = " << procID() << endl;
  unpackReceivesToMe(a_op); // nullOp in uniprocessor mode
#endif
}
#ifndef PR_MPI
template<typename T, unsigned int C, MemType MEM, unsigned int CTR>
void LevelBoxData<T, C, MEM, CTR>::writeSendDataFromMeIntoBuffers(const LDOP& a_op) const
{
}
template<typename T, unsigned int C, MemType MEM, unsigned int CTR>
void LevelBoxData<T, C, MEM, CTR>::postSendsFromMe() const
{
}
template<typename T, unsigned int C, MemType MEM, unsigned int CTR>
void LevelBoxData<T, C, MEM, CTR>::postReceivesToMe() const
{
}
template<typename T, unsigned int C, MemType MEM, unsigned int CTR>
void LevelBoxData<T, C, MEM, CTR>::unpackReceivesToMe(const LDOP& a_op) const
{
}
#else
//========================================================================
//
// data structures used by makeItSo when we have some
// data that needs to be moved (ie. there are entries
// in the 'FROM' or 'TO' CopyIterators)
//
template<typename T, unsigned int C, MemType MEM, unsigned int CTR>
void 
LevelBoxData<T, C, MEM, CTR>::completePendingSends() const
{
  PR_TIME("completePendingSends");
  //pout() << "numSends " << this->m_buff->numSends << endl;
  if (this->m_buff->numSends > 0)
    {   
      PR_TIME("MPI_Waitall");
      m_buff->m_sendStatus.resize(this->m_buff->numSends);
      //pout() << "this->m_buff->numSends " << this->m_buff->numSends << ", proc = " << procID() << endl;
      //pout() << "&(m_buff->m_sendRequests[0]) " << &(m_buff->m_sendRequests[0]) <<
      // " , " <<  m_buff->m_sendRequests.size() << endl;
        //pout() << "&(m_buff->m_sendStatus[0]) " << &(m_buff->m_sendStatus[0]) <<
      // " , " << m_buff->m_sendStatus.size() << endl; 
      int result = MPI_Waitall(this->m_buff->numSends, &(m_buff->m_sendRequests[0]), &(m_buff->m_sendStatus[0]));
      //pout() << "finished MPI_Waitall" << endl;
      if (result != MPI_SUCCESS)
        {
          cout << "result - failed " << result << endl;
          //hell if I know what to do about failed messaging here
        }
    }
  this->m_buff->numSends = 0;
}
template<typename T, unsigned int C, MemType MEM, unsigned int CTR>
void // Getting rid of component dependencies.
LevelBoxData<T, C, MEM, CTR>::allocateBuffers(
                     const Copier&   a_copier,
                     const LDOP& a_op) const
{
  PR_TIME("MPI_allocateBuffers");
  m_buff = &(((Copier&)a_copier).m_buffers);
  
  if (m_buff->isDefined()) return;
  
  m_buff->m_fromMe.resize(0);
  m_buff->m_toMe.resize(0);
  size_t sendBufferSize = 0;
  size_t recBufferSize  = 0;

  BoxData<T,C> dummy;
  for (CopyIterator it(a_copier, CopyIterator::FROM); it.ok(); ++it)
    {
      const MotionItem& item = it();
      CopierBuffer::bufEntry b;
      b.item = &item;
      b.size = a_op.size(item.fromRegion);
      sendBufferSize+=b.size;
      b.procID = item.procID;
      m_buff->m_fromMe.push_back(b);
    }
  sort(m_buff->m_fromMe.begin(), m_buff->m_fromMe.end());
  for (CopyIterator it(a_copier, CopyIterator::TO); it.ok(); ++it)
    {
      const MotionItem& item = it();
      CopierBuffer::bufEntry b;
      b.item = &item;
      b.size = a_op.size(item.fromRegion);
      recBufferSize+=b.size;
      b.procID = item.procID;
      m_buff->m_toMe.push_back(b);
    }
  sort(m_buff->m_toMe.begin(), m_buff->m_toMe.end());

  // allocate send and receive buffer space.

  if (sendBufferSize > m_buff->m_sendcapacity)
    {
      free((m_buff->m_sendbuffer));
      if (s_verbosity > 0) pout()<<"malloc send buffer "<<sendBufferSize<<std::endl;
      (m_buff->m_sendbuffer) = malloc(sendBufferSize);
      if ((m_buff->m_sendbuffer) == NULL)
        {
          MayDay<void>::Error("Out of memory in BoxLayoutData::allocatebuffers");
        }
      m_buff->m_sendcapacity = sendBufferSize;
    }

  if (recBufferSize > m_buff->m_reccapacity)
    {
      free(m_buff->m_recbuffer);
      if (s_verbosity > 0) pout()<<"malloc receive buffer "<<recBufferSize<<std::endl;
      m_buff->m_recbuffer = malloc(recBufferSize);
      if (m_buff->m_recbuffer == NULL)
        {
          MayDay<void>::Error("Out of memory in BoxLayoutData::allocatebuffers");
        }
      m_buff->m_reccapacity = recBufferSize;
    }

  char* nextFree = (char*)(m_buff->m_sendbuffer);
  if (m_buff->m_fromMe.size() > 0)
    {
      for (unsigned int i=0; i<m_buff->m_fromMe.size(); ++i)
        {
          m_buff->m_fromMe[i].bufPtr = nextFree;
          nextFree += m_buff->m_fromMe[i].size;
        }
    }

  nextFree = (char*)m_buff->m_recbuffer;
  if (m_buff->m_toMe.size() > 0)
    {
      for (unsigned int i=0; i<m_buff->m_toMe.size(); ++i)
        {
          m_buff->m_toMe[i].bufPtr = nextFree;
          nextFree += m_buff->m_toMe[i].size;
        }
    }
  // since fromMe and toMe are sorted based on procID, messages can now be grouped
  // together on a per-processor basis.
}
template<typename T, unsigned int C, MemType MEM, unsigned int CTR> 
void 
LevelBoxData<T, C, MEM, CTR>::writeSendDataFromMeIntoBuffers(const LDOP& a_op) const
{
  PR_TIME("write Data to buffers");
  int isize = m_buff->m_fromMe.size();
  for (unsigned int i=0; i< isize; ++i)
    {
      const CopierBuffer::bufEntry& entry = m_buff->m_fromMe[i];
      a_op.linearOut(entry.bufPtr,entry.item->fromRegion,entry.item->fromIndex);
    }
}
template<typename T, unsigned int C, MemType MEM, unsigned int CTR>
void 
LevelBoxData<T, C, MEM, CTR>::postSendsFromMe() const
{
  PR_TIME("post_Sends");
  // now we get the magic of message coalescence
  // fromMe has already been sorted in the allocateBuffers() step.
  // PC: Looks like we are just compressing out the buffers by setting the size to zero.
  // This is called after allocateBuffers.
  this->m_buff->numSends = m_buff->m_fromMe.size();
  if (this->m_buff->numSends > 1)
    {
      for (unsigned int i=m_buff->m_fromMe.size()-1; i>0; --i)
        {
          if (m_buff->m_fromMe[i].procID == m_buff->m_fromMe[i-1].procID)
            {
              // Decrementing the number of sends as we coalesce.
              this->m_buff->numSends--;
              m_buff->m_fromMe[i-1].size = m_buff->m_fromMe[i-1].size + m_buff->m_fromMe[i].size;
              m_buff->m_fromMe[i].size = 0;
            }
        }
    }
  m_buff->m_sendRequests.resize(this->m_buff->numSends);
  std::list<MPI_Request> extraRequests;
    
  unsigned int next=0;
  long long maxSize = 0;
  for (int i=0; i<this->m_buff->numSends; ++i)
    {
      const CopierBuffer::bufEntry& entry = m_buff->m_fromMe[next];
      char*  buffer = (char*)entry.bufPtr;
      std::size_t bsize = entry.size;
      int idtag=0;
      while (bsize > PR_MAX_MPI_MESSAGE_SIZE)
        {
          extraRequests.push_back(MPI_Request());
          {
            PR_TIME("MPI_Isend");
            //pout() << "postSendsFromMe -MPI_Isend extra" << endl;
            MPI_Isend(buffer, PR_MAX_MPI_MESSAGE_SIZE, MPI_BYTE, entry.procID,
                      idtag, Proto_MPI<void>::comm, &(extraRequests.back()));
          }
          maxSize = PR_MAX_MPI_MESSAGE_SIZE;
          bsize -= PR_MAX_MPI_MESSAGE_SIZE;
          buffer+=PR_MAX_MPI_MESSAGE_SIZE;
          idtag++;
        }
      {
        PR_TIME("MPI_Isend");
        //pout() << "postSendsFromMe - MPI_Isend main" << endl;
        MPI_Isend(buffer, bsize, MPI_BYTE, entry.procID,
                  idtag, Proto_MPI<void>::comm, &(m_buff->m_sendRequests[i]));
      }
      maxSize = std::max<long long>(bsize, maxSize);
      ++next;
      while (next < m_buff->m_fromMe.size() && m_buff->m_fromMe[next].size == 0) ++next;
    }
  for (std::list<MPI_Request>::iterator it = extraRequests.begin(); it != extraRequests.end(); ++it)
    {
      m_buff->m_sendRequests.push_back(*it);
    }
  this->m_buff->numSends = m_buff->m_sendRequests.size();
    
  //PR_MaxMPISendSize = std::max<long long>(PR_MaxMPISendSize, maxSize);
}
template<typename T, unsigned int C, MemType MEM, unsigned int CTR>
void 
LevelBoxData<T, C, MEM, CTR>::postReceivesToMe() const
{
  PR_TIME("post_Receives");
  this->m_buff->numReceives = m_buff->m_toMe.size();
  
  if (this->m_buff->numReceives > 1)
    {
      for (unsigned int i=m_buff->m_toMe.size()-1; i>0; --i)
        {
          if (m_buff->m_toMe[i].procID == m_buff->m_toMe[i-1].procID)
            {
              this->m_buff->numReceives--;
              m_buff->m_toMe[i-1].size += m_buff->m_toMe[i].size;
              m_buff->m_toMe[i].size = 0;
            }
          
        }
    }
  m_buff->m_receiveRequests.resize(this->m_buff->numReceives);
  std::list<MPI_Request> extraRequests;
  unsigned int next=0;
  long long maxSize = 0;
  for (int i=0; i<this->m_buff->numReceives; ++i)
    {
      const CopierBuffer::bufEntry& entry = m_buff->m_toMe[next];
      char*  buffer = (char*)entry.bufPtr;
      size_t bsize = entry.size;
      int idtag=0;
      while (bsize > PR_MAX_MPI_MESSAGE_SIZE)
        {
          extraRequests.push_back(MPI_Request());
          {
            PR_TIME("MPI_Irecv");
            // pout() << "PostReceivesToMe - MPI_Irecv - extra" << endl;
            MPI_Irecv(buffer,PR_MAX_MPI_MESSAGE_SIZE, MPI_BYTE, entry.procID,
                      idtag, Proto_MPI<void>::comm, &(extraRequests.back()));
          }
          maxSize = PR_MAX_MPI_MESSAGE_SIZE;
          bsize -= PR_MAX_MPI_MESSAGE_SIZE;
          buffer+=PR_MAX_MPI_MESSAGE_SIZE;
          idtag++;
        }
      {
        PR_TIME("MPI_Irecv");
        //pout() << "PostReceivesToMe - MPI_Irecv - main" << endl;
        MPI_Irecv(buffer, bsize, MPI_BYTE, entry.procID,
                  idtag, Proto_MPI<void>::comm, &(m_buff->m_receiveRequests[i]));
      }
      ++next;
      maxSize = std::max<long long>(bsize, maxSize);
      while (next < m_buff->m_toMe.size() && m_buff->m_toMe[next].size == 0) ++next;
    }
  for (std::list<MPI_Request>::iterator it = extraRequests.begin(); it != extraRequests.end(); ++it)
    {
      m_buff->m_receiveRequests.push_back(*it);
    }
  this->m_buff->numReceives = m_buff->m_receiveRequests.size();

  // PR_MaxMPIRecvSize = std::max<long long>(PR_MaxMPIRecvSize, maxSize);
  //pout()<<"maxSize="<<maxSize<<" posted "<<this->m_buff->numReceives<<" receives\n";
}
template<typename T, unsigned int C, MemType MEM, unsigned int CTR>
void 
LevelBoxData<T, C, MEM, CTR>::unpackReceivesToMe(const LDOP& a_op) const
{
  PR_TIME("unpack_messages");
  // pout() << "numReceives " << this->m_buff->numReceives << endl;
  if (this->m_buff->numReceives > 0)
    {
      m_buff->m_receiveStatus.resize(this->m_buff->numReceives);
      int result;
      {
        PR_TIME("MPI_Waitall");
        //pout() << "this->m_buff->numReceives " << this->m_buff->numReceives << endl;
        //pout() << "size (m_buff->m_receiveRequests[0]) " << &(m_buff->m_receiveRequests[0]) <<
        //  " , " <<  m_buff->m_receiveRequests.size() << endl;
        //pout() << "&(m_buff->m_receiveStatus[0]) " << &(m_buff->m_receiveStatus[0]) <<
        //  " , " << m_buff->m_receiveStatus.size() << endl; 
        result = MPI_Waitall(this->m_buff->numReceives, &(m_buff->m_receiveRequests[0]),
                             &(m_buff->m_receiveStatus[0]));
        //pout() << "unpackReceives - post Waitall " << this->m_buff->numReceives << endl;
      }
      if (result != MPI_SUCCESS)
        {
          //hell if I know what to do about failed messaging here
          //maybe a mayday::warning?
        }

      int isize = m_buff->m_toMe.size();
      for (unsigned int i=0; i< isize; ++i)
        {
          const CopierBuffer::bufEntry& entry = m_buff->m_toMe[i];
          a_op.linearIn(entry.bufPtr,entry.item->toRegion,entry.item->toIndex);
        }
    }
  this->m_buff->numReceives = 0;
}
#endif
*/
#if 0
void WriteData(LevelBoxData<double>& a_phi,int a_iter,double a_dx,std::string a_str)
{
  DisjointBoxLayout bl = a_phi.getDBL();
  Box bx = bl.getDomain();
  BoxData<double> phiOut(bx);
  for (int i = 0; i < bl.size();i++)
    {
      a_phi[i].copyTo(phiOut,bl[i]);
    }
  WriteData(phiOut,a_iter,a_dx,a_str.c_str(),a_str.c_str());
}
#endif
