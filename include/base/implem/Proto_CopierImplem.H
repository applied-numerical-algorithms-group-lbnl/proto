//========================================================================
// COPIER PUBLIC API

template<class OP, typename P_SRC, typename P_DST, MemType SRC_MEM, MemType DST_MEM>
Copier<OP, P_SRC, P_DST, SRC_MEM, DST_MEM>::Copier()
{
    m_sendBuffer = nullptr;
    m_recvBuffer  = nullptr;
    m_sendCapacity = 0;
    m_recvCapacity = 0;
#ifdef PR_MPI
    m_fromMe.clear();
    m_toMe.clear();
    m_numSends = 0;
    m_numRecvs = 0;
#endif
    m_isDefined = false;
}

template<class OP, typename P_SRC, typename P_DST, MemType SRC_MEM, MemType DST_MEM>
Copier<OP, P_SRC, P_DST, SRC_MEM, DST_MEM>::Copier(OP a_op)
{
    define(a_op);
}

template<class OP, typename P_SRC, typename P_DST, MemType SRC_MEM, MemType DST_MEM>
Copier<OP, P_SRC, P_DST, SRC_MEM, DST_MEM>::~Copier()
{
    clearBuffers();
}

template<class OP, typename P_SRC, typename P_DST, MemType SRC_MEM, MemType DST_MEM>
void
Copier<OP, P_SRC, P_DST, SRC_MEM, DST_MEM>::define(OP a_op)
{
    clearBuffers();
    clearMotionPlan();
    m_op = a_op;
    this->buildMotionPlans(m_op); 
    m_isDefined = true;
}

template<class OP, typename P_SRC, typename P_DST, MemType SRC_MEM, MemType DST_MEM>
void
Copier<OP, P_SRC, P_DST, SRC_MEM, DST_MEM>::clearBuffers()
{
    if (m_sendBuffer != nullptr) { proto_free<SRC_MEM>(m_sendBuffer); }
    if (m_recvBuffer != nullptr) { proto_free<DST_MEM>(m_recvBuffer); }
    m_sendBuffer = nullptr;
    m_recvBuffer  = nullptr;
    m_sendCapacity = 0;
    m_recvCapacity = 0;
#ifdef PR_MPI
    m_fromMe.clear();
    m_toMe.clear();
    m_numSends = 0;
    m_numRecvs = 0;
#endif
}

template<class OP, typename P_SRC, typename P_DST, MemType SRC_MEM, MemType DST_MEM>
void
Copier<OP, P_SRC, P_DST, SRC_MEM, DST_MEM>::clear()
{
    clearBuffers();
}

template<class OP, typename P_SRC, typename P_DST, MemType SRC_MEM, MemType DST_MEM>
void
Copier<OP, P_SRC, P_DST, SRC_MEM, DST_MEM>::clearMotionPlan()
{
    m_localMotionPlan.clear();
    m_fromMotionPlan.clear();
    m_toMotionPlan.clear();

    m_isDefined = false;
}

template<class OP, typename P_SRC, typename P_DST, MemType SRC_MEM, MemType DST_MEM>
void
Copier<OP, P_SRC, P_DST, SRC_MEM, DST_MEM>::sort()
{
    std::sort(m_fromMotionPlan.begin(), m_fromMotionPlan.end());
    std::sort(m_toMotionPlan.begin(), m_toMotionPlan.end());
}

template<class OP, typename P_SRC, typename P_DST, MemType SRC_MEM, MemType DST_MEM>
std::vector<MotionItem<P_SRC,P_DST>>&
Copier<OP, P_SRC, P_DST, SRC_MEM, DST_MEM>::motionPlan(MotionType a_type)
{
    switch (a_type)
    {
        case LOCAL: return m_localMotionPlan;
        case FROM:  return m_fromMotionPlan;
        case TO:    return m_toMotionPlan;
    }
}


template<class OP, typename P_SRC, typename P_DST, MemType SRC_MEM, MemType DST_MEM>
bool
Copier<OP, P_SRC, P_DST, SRC_MEM, DST_MEM>::operator==(const Copier<OP, P_SRC, P_DST, SRC_MEM, DST_MEM>& a_rhs) const
{
    if (m_localMotionPlan != a_rhs.m_localMotionPlan) {return false;}
    if (m_fromMotionPlan != a_rhs.m_fromMotionPlan) {return false;}
    if (m_toMotionPlan != a_rhs.m_toMotionPlan) {return false;}
}

template<class OP, typename P_SRC, typename P_DST, MemType SRC_MEM, MemType DST_MEM>
void
Copier<OP, P_SRC, P_DST, SRC_MEM, DST_MEM>::execute()
{
    barrier();
    // execute() and makeItSo() are aliases of eachother.
    // "makeItSo" is maintained for the sake of posterity
    // and clarity to the development team who are used to
    // the semantics of this function -CLG
    makeItSo();
}

template<class OP, typename P_SRC, typename P_DST, MemType SRC_MEM, MemType DST_MEM>
CopierIterator<P_SRC,P_DST>
Copier<OP, P_SRC, P_DST, SRC_MEM, DST_MEM>::begin(MotionType a_type) const
{
    switch (a_type)
    {
        case FROM:
            return CopierIterator<P_SRC,P_DST>(&m_fromMotionPlan);
        case TO:
            return CopierIterator<P_SRC,P_DST>(&m_toMotionPlan);
        case LOCAL:
            return CopierIterator<P_SRC,P_DST>(&m_localMotionPlan);
    }
    MayDay<void>::Abort("Copier::begin | Error: invalid MotionType");
    return 0;
}
//========================================================================
// COPIER PRIVATE FUNCTIONS

template<class OP, typename P_SRC, typename P_DST, MemType SRC_MEM, MemType DST_MEM>
void
Copier<OP, P_SRC, P_DST, SRC_MEM, DST_MEM>::makeItSo()
{
    makeItSoBegin();
    makeItSoLocal();
    makeItSoEnd();
}

template<class OP, typename P_SRC, typename P_DST, MemType SRC_MEM, MemType DST_MEM>
void
Copier<OP, P_SRC, P_DST, SRC_MEM, DST_MEM>::makeItSoBegin()
{
#ifdef PR_MPI
    PR_TIME("Copier::makeItSoBegin");
    allocateBuffers();
    writeToSendBuffers();
    m_numRecvs = m_toMe.size();
    if (m_numRecvs > 0)
    {
        postRecvs(); // non-blocking
    }
    m_numSends = m_fromMe.size();
    if (m_numSends > 0)
    {
        postSends(); // non-blocking
    }
#endif
}

template<class OP, typename P_SRC, typename P_DST, MemType SRC_MEM, MemType DST_MEM>
void
Copier<OP, P_SRC, P_DST, SRC_MEM, DST_MEM>::makeItSoLocal()
{
    PR_TIME("Copier::makeItSoLocal");
    for (auto iter = begin(LOCAL); iter.ok(); ++iter)
    {
        const auto& item = *iter;
        //m_op.localCopy(item.fromRegion, item.fromIndex, item.toRegion, item.toIndex);
        m_op.localCopy(item);
    }
}

template<class OP, typename P_SRC, typename P_DST, MemType SRC_MEM, MemType DST_MEM>
void
Copier<OP, P_SRC, P_DST, SRC_MEM, DST_MEM>::makeItSoEnd()
{
#ifdef PR_MPI
    PR_TIME("Copier::makeItSoEnd");
   completePendingSends();
   readFromRecvBuffers(); 
#endif
}

template<class OP, typename P_SRC, typename P_DST, MemType SRC_MEM, MemType DST_MEM>
void
Copier<OP, P_SRC, P_DST, SRC_MEM, DST_MEM>::completePendingSends() const
{
#ifdef PR_MPI
    PR_TIME("Copier::completePendingSends");
    if (m_numSends > 0)
    {
        m_sendStatus.resize(m_numSends);
    }
    int result = MPI_Waitall(m_numSends, &m_sendRequests[0], &m_sendStatus[0]);
    if (result != MPI_SUCCESS)
    {
        // TODO: Figure out what to do here. 
    }
    m_numSends = 0;
#endif
}

template<class OP, typename P_SRC, typename P_DST, MemType SRC_MEM, MemType DST_MEM>
void
Copier<OP, P_SRC, P_DST, SRC_MEM, DST_MEM>::allocateBuffers()
{
#ifdef PR_MPI
    // clear buffers
    m_fromMe.resize(0);
    m_toMe.resize(0);
    size_t sendBufferSize = 0;
    size_t recvBufferSize = 0;

    // read motion plans
    for (auto iter = begin(FROM); iter.ok(); ++iter)
    {
        auto& item = *iter;
        BufferEntry<P_SRC,P_DST> b;
        b.item = &item;
        b.size = m_op.linearSize(item.fromRegion, item.fromIndex);
        sendBufferSize += b.size;
        b.procID = item.procID;
        m_fromMe.push_back(b);
    }
    std::sort(m_fromMe.begin(), m_fromMe.end());
    for (auto iter = begin(TO); iter.ok(); ++iter)
    {
        auto& item = *iter;
        BufferEntry<P_SRC,P_DST> b;
        b.item = &item;
        b.size = m_op.linearSize(item.fromRegion, item.fromIndex);
        recvBufferSize += b.size;
        b.procID = item.procID;
        m_toMe.push_back(b);
    }
    std::sort(m_toMe.begin(), m_toMe.end());
    // allocate send and recv buffers if needed
    if (sendBufferSize > m_sendCapacity)
    {
        if (m_sendCapacity > 0) { proto_free<SRC_MEM>(m_sendBuffer); }
        m_sendBuffer = proto_malloc<SRC_MEM>(sendBufferSize);
        if (m_sendBuffer == NULL)
        {
            MayDay<void>::Error("Copier::allocateBuffers | Error: Out of memory.");
        }
        m_sendCapacity = sendBufferSize;
    }
    if (recvBufferSize > m_recvCapacity)
    {
        if (m_recvCapacity > 0) { proto_free<DST_MEM>(m_recvBuffer); }
        m_recvBuffer = proto_malloc<DST_MEM>(recvBufferSize);
        if (m_recvBuffer == NULL)
        {
            MayDay<void>::Error("Copier::allocateBuffers | Error: Out of memory.");
        }
        m_recvCapacity = recvBufferSize;
    }

    char* nextFree = (char*)m_sendBuffer;
    if (this->m_fromMe.size() > 0)
    {
        for (unsigned int ii = 0; ii < m_fromMe.size(); ii++)
        {
            m_fromMe[ii].buffer = nextFree;
            nextFree += m_fromMe[ii].size;
        }
    }
    
    nextFree = (char*)m_recvBuffer;
    if (this->m_toMe.size() > 0)
    {
        for (unsigned int ii = 0; ii < m_toMe.size(); ii++)
        {
            m_toMe[ii].buffer = nextFree;
            nextFree += m_toMe[ii].size;
        }
    }
#endif
}

template<class OP, typename P_SRC, typename P_DST, MemType SRC_MEM, MemType DST_MEM>
void
Copier<OP, P_SRC, P_DST, SRC_MEM, DST_MEM>::writeToSendBuffers()
{
#ifdef PR_MPI
    PR_TIME("Copier::writeToSendBuffers");
    for (unsigned int ii = 0; ii < m_fromMe.size(); ii++)
    {
        auto item = m_fromMe[ii];
        m_op.linearOut(item.buffer, *item.item);
    }
    barrier();
#endif
}

template<class OP, typename P_SRC, typename P_DST, MemType SRC_MEM, MemType DST_MEM>
void
Copier<OP, P_SRC, P_DST, SRC_MEM, DST_MEM>::postSends() const
{
#ifdef PR_MPI
    PR_TIME("Copier::postSends");
    // Coalesce messages
    m_numSends = m_fromMe.size();
    if (m_numSends > 1)
    {
        for (unsigned int ii = m_fromMe.size()-1; ii > 0; ii--)
        {
            if (m_fromMe[ii].procID == m_fromMe[ii-1].procID)
            {
                m_numSends--;
                m_fromMe[ii-1].size += m_fromMe[ii].size;
                m_fromMe[ii].size = 0;
            }
        }
    }
    
    // post the send requests
    m_sendRequests.resize(m_numSends);
    std::list<MPI_Request> extraRequests;
    unsigned int next = 0;
    //long long maxSize = 0; //I don't think maxSize is used for anything -CLG
    for (int ii = 0; ii < m_numSends; ii++)
    {
        const auto& entry = m_fromMe[next];
        char* buffer = (char*)entry.buffer;
        std::size_t bufferSize = entry.size;
        int idtag = 0;
        // deal with possible message overflow
        while (bufferSize > PR_MAX_MPI_MESSAGE_SIZE)
        {
            extraRequests.push_back(MPI_Request());
            {
                PR_TIME("MPI_Isend");
                MPI_Isend(buffer, PR_MAX_MPI_MESSAGE_SIZE, MPI_BYTE,
                    entry.procID, idtag, Proto_MPI<void>::comm,
                    &(extraRequests.back()));
            }
            //maxSize = PR_MAX_MPI_MESSAGE_SIZE;
            bufferSize -= PR_MAX_MPI_MESSAGE_SIZE;
            buffer += PR_MAX_MPI_MESSAGE_SIZE;
            idtag++;
        }
        {
            PR_TIME("MPI_Isend");
            MPI_Isend(buffer, bufferSize, MPI_BYTE, entry.procID,
                idtag, Proto_MPI<void>::comm, &(m_sendRequests[ii]));
        }
        //maxSize = std::max<long long>(bufferSize, maxSize);
        next++;
        // skip all the size==0 messages from coalescence
        while (next < m_fromMe.size() && m_fromMe[next].size == 0) { next++; }
    }
    for (auto iter = extraRequests.begin(); iter != extraRequests.end(); ++iter)
    {
        m_sendRequests.push_back(*iter);
    }
    m_numSends = m_sendRequests.size();
#endif
}

template<class OP, typename P_SRC, typename P_DST, MemType SRC_MEM, MemType DST_MEM>
void
Copier<OP, P_SRC, P_DST, SRC_MEM, DST_MEM>::readFromRecvBuffers()
{
#ifdef PR_MPI
    PR_TIME("Copier::readFromRecvBuffers");
    if (m_numRecvs > 0)
    {
        m_recvStatus.resize(m_numRecvs);
        int result;
        {
            PR_TIME("MPI_Waitall");
            result = MPI_Waitall(m_numRecvs, &(m_recvRequests[0]),
                &(m_recvStatus[0]));
        }
        if (result != MPI_SUCCESS)
        {
            // TODO: figure out what to do here
        }
        
        for (unsigned int ii = 0; ii < m_toMe.size(); ii++)
        {
            const auto& entry = m_toMe[ii];
            
            m_op.linearIn(entry.buffer, *entry.item);
        }
    }
    m_numRecvs = 0;
#endif
}

template<class OP, typename P_SRC, typename P_DST, MemType SRC_MEM, MemType DST_MEM>
void
Copier<OP, P_SRC, P_DST, SRC_MEM, DST_MEM>::postRecvs() const
{
#ifdef PR_MPI
    PR_TIME("Copier::postRecvs");
    m_numRecvs = m_toMe.size();
    // Do message coalescence
    if (m_numRecvs > 1)
    {
        for (unsigned int ii = m_toMe.size()-1; ii > 0; ii--)
        {
            if (m_toMe[ii].procID == m_toMe[ii-1].procID)
            {
                m_numRecvs--;
                m_toMe[ii-1].size += m_toMe[ii].size;
                m_toMe[ii].size = 0;
            }
        }
    }
    m_recvRequests.resize(m_numRecvs);
    std::list<MPI_Request> extraRequests;
    unsigned int next = 0;
    //long long maxSize = 0; // I don't think this gets used -CLG
    for (int ii = 0; ii < m_numRecvs; ii++)
    {
        const auto& entry = m_toMe[next];
        char* buffer = (char*)entry.buffer;
        size_t bufferSize = entry.size;
        int idtag = 0;
        while (bufferSize > PR_MAX_MPI_MESSAGE_SIZE)
        {
            extraRequests.push_back(MPI_Request());
            {
                PR_TIME("MPI_Irecv");
                MPI_Irecv(buffer, PR_MAX_MPI_MESSAGE_SIZE, MPI_BYTE, entry.procID,
                        idtag, Proto_MPI<void>::comm, &(extraRequests.back()));
            }
            bufferSize -= PR_MAX_MPI_MESSAGE_SIZE;
            buffer += PR_MAX_MPI_MESSAGE_SIZE;
            idtag++;
        }
        {
            PR_TIME("MPI_Irecv");
            MPI_Irecv(buffer, bufferSize, MPI_BYTE, entry.procID,
                idtag, Proto_MPI<void>::comm, &(m_recvRequests[ii]));
        }
        next++;
        while(next < m_toMe.size() && m_toMe[next].size == 0) { ++next; }
    }
    for (auto iter = extraRequests.begin(); iter != extraRequests.end(); ++iter)
    {
        m_recvRequests.push_back(*iter);
    }
    m_numRecvs = m_recvRequests.size();
#endif
}

template<class OP, typename P_SRC, typename P_DST, MemType SRC_MEM, MemType DST_MEM>
void
Copier<OP, P_SRC, P_DST, SRC_MEM, DST_MEM>::print() const
{
    using Proto::pr_out;
    pr_out() << std::endl << "-----------------------------------------------------------------------" << std::endl;
    pr_out() << "Local Copies: " << std::endl;
    for (auto iter = m_localMotionPlan.begin(); iter != m_localMotionPlan.end(); ++iter)
    {
        auto item = *iter;
        pr_out() << item.fromRegion << " @ " << item.fromIndex.global();
        pr_out() << " ==> " << item.toRegion << " @ " << item.toIndex.global() << std::endl;
    }
    pr_out() << "Outgoing Data: " << std::endl;
    for (auto iter = m_fromMotionPlan.begin(); iter != m_fromMotionPlan.end(); ++iter)
    {
        auto item = *iter;
        pr_out() << item.fromRegion << " @ " << item.fromIndex.global();
        pr_out() << " ==> " << item.toRegion << " @ " << item.toIndex.global() << std::endl;
    }
    pr_out() << "Incoming Data: " << std::endl;
    for (auto iter = m_toMotionPlan.begin(); iter != m_toMotionPlan.end(); ++iter)
    {
        auto item = *iter;
        pr_out() << item.fromRegion << " @ " << item.fromIndex.global();
        pr_out() << " ==> " << item.toRegion << " @ " << item.toIndex.global() << std::endl;
    }
}
//========================================================================
// MOTIONITEM, BUFFERENTRY, and COPYITERATOR FUNCTIONS

template<typename P_SRC, typename P_DST>
MotionItem<P_SRC,P_DST>::MotionItem(const DataIndex<P_SRC>& a_fromIndex,
        const DataIndex<P_DST>& a_toIndex,
        const Box& a_fromRegion,
        const Box& a_toRegion)
{
    PROTO_ASSERT(a_fromRegion.size() == a_toRegion.size(),
            "MotionItem::Constructor | toRegion and fromRegion should be the same size.");
    fromIndex  = a_fromIndex;
    toIndex    = a_toIndex;
    fromRegion = a_fromRegion;
    toRegion   = a_toRegion;
    procID     = Proto::procID();
    // rotation defaults to identity
}

template<typename P_SRC, typename P_DST>
bool
MotionItem<P_SRC,P_DST>::operator==(const MotionItem<P_SRC,P_DST>& a_rhs) const
{
    return (fromIndex  == a_rhs.fromIndex) &&
           (toIndex    == a_rhs.toIndex) &&
           (fromRegion == a_rhs.fromRegion) &&
           (toRegion   == a_rhs.toRegion) &&
           (procID     == a_rhs.procID) &&
           (rotation   == a_rhs.rotation);
}

template<typename P_SRC, typename P_DST>
bool
MotionItem<P_SRC,P_DST>::operator<(const MotionItem<P_SRC,P_DST>& a_rhs) const
{
    if (procID != a_rhs.procID)
    {
        return (procID < a_rhs.procID);
    } else if (toIndex.global() != a_rhs.toIndex.global())
    {
        return (toIndex.global() < a_rhs.toIndex.global());
    } else if (fromIndex.global() != a_rhs.fromIndex.global())
    {
        return (fromIndex.global() < a_rhs.fromIndex.global());
    } else if (toRegion != a_rhs.toRegion)
    {
        return (toRegion < a_rhs.toRegion);
    } else
    {
        return (fromRegion < a_rhs.fromRegion);
    }
    // rotation is ignored
}

template<typename P_SRC, typename P_DST>
bool
BufferEntry<P_SRC,P_DST>::operator<(const BufferEntry<P_SRC,P_DST>& a_rhs) const
{
    return ((*item) < *(a_rhs.item));
}

template<typename P_SRC, typename P_DST>
CopierIterator<P_SRC,P_DST>::CopierIterator(const std::vector<MotionItem<P_SRC,P_DST>>* a_motionPlan)
{
    m_current = 0;
    m_motionPlan = a_motionPlan;
}

template<typename P_SRC, typename P_DST>
const MotionItem<P_SRC, P_DST>&
CopierIterator<P_SRC, P_DST>::operator*() const
{
    PROTO_ASSERT(m_current < size(),
        "CopierIterator::operator() | Error: Index out of bounds: %i >= %lu.",
        m_current, size());
    return m_motionPlan->operator[](m_current);
}

template<typename P_SRC, typename P_DST>
const MotionItem<P_SRC, P_DST>&
CopierIterator<P_SRC, P_DST>::operator()() const
{
    return *(*this);
}

template<typename P_SRC, typename P_DST>
const MotionItem<P_SRC, P_DST>&
CopierIterator<P_SRC, P_DST>::operator[](size_t a_index) const
{
    PROTO_ASSERT(a_index < size(),
        "CopierIterator::operator() | Error: Index out of bounds: %lu >= %lu.",
        a_index, size());
    return m_motionPlan->operator[](a_index);
}

template<typename P_SRC, typename P_DST>
void
CopierIterator<P_SRC, P_DST>::operator++()
{
    ++m_current;
}

template<typename P_SRC, typename P_DST>
bool
CopierIterator<P_SRC, P_DST>::ok()
{
    return (m_current < size());
}

template<typename P_SRC, typename P_DST>
void
CopierIterator<P_SRC, P_DST>::reset()
{
    m_current = 0;
}

template<typename P_SRC, typename P_DST>
size_t
CopierIterator<P_SRC, P_DST>::size() const
{
    return m_motionPlan->size();
}


