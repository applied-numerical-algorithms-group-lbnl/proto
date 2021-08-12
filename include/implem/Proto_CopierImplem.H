//========================================================================
// COPIER PUBLIC API

template<class OP>
Copier<OP>::Copier()
{
}

template<class OP>
Copier<OP>::Copier(OP a_op)
{
    define(a_op);
}

template<class OP>
Copier<OP>::~Copier()
{
    if (m_isDefined)
    {
        clear();
    }
}

template<class OP>
void
Copier<OP>::define(OP a_op)
{
    m_op = a_op;
    m_numSends = 0;
    m_numRecvs = 0;
    m_sendCapacity = 0;
    m_recvCapacity = 0;
   // m_sendBuffer = NULL;
   //m_recvBuffer = NULL;
    this->buildMotionPlans(m_op); 
    m_isDefined = true;
}

template<class OP>
void
Copier<OP>::clear()
{
    //if (m_sendBuffer != NULL) { free(m_sendBuffer); }
    //if (m_recvBuffer != NULL) { free(m_recvBuffer); }
    //m_sendBuffer = NULL;
    //m_recvBuffer  = NULL;
    m_sendCapacity = 0;
    m_recvCapacity = 0;
   
    // TODO: clear the motion plan arrays 
    
    m_isDefined = false;
}

template<class OP>
void
Copier<OP>::sort()
{
    std::vector<MotionItem*>& vfrom  = m_fromMotionPlan;
    std::sort(vfrom.begin(), vfrom.end(), MotionItemSorter());
    std::vector<MotionItem*>& vto = m_toMotionPlan;
    std::sort(vto.begin(), vto.end(), MotionItemSorter());

}

template<class OP>
Copier<OP>&
Copier<OP>::operator=(const Copier<OP>& a_rhs)
{
}

template<class OP>
bool
Copier<OP>::operator==(const Copier<OP>& a_rhs) const
{
    // Because we store *pointers* to MotionItem, we can't just compare the
    // std::vectors that define the motion plans. I see no reason to do this
    // since MotionItem is always the same size. -CLG
    if (m_localMotionPlan.size() != a_rhs.m_localMotionPlan.size()) { return false;}
    if (m_fromMotionPlan.size()  != a_rhs.m_fromMotionPlan.size())  { return false;}
    if (m_toMotionPlan.size()    != a_rhs.m_toMotionPlan.size())    { return false;}
    for (int ii = 0; ii < m_localMotionPlan.size(); ii++)
    {
        if (!((*m_localMotionPlan[ii]) == (*a_rhs.m_localMotionPlan[ii])))
        {
            return false;
        }
    }
    for (int ii = 0; ii < m_fromMotionPlan.size(); ii++)
    {
        if (!((*m_fromMotionPlan[ii]) == (*a_rhs.m_fromMotionPlan[ii])))
        {
            return false;
        }
    }
    for (int ii = 0; ii < m_toMotionPlan.size(); ii++)
    {
        if (!((*m_toMotionPlan[ii]) == (*a_rhs.m_toMotionPlan[ii])))
        {
            return false;
        }
    }
}

template<class OP>
void
Copier<OP>::reverse()
{
}

template<class OP>
void
Copier<OP>::execute()
{
    // execute() and makeItSo() are aliases of eachother.
    // "makeItSo" is maintained for the sake of posterity
    // and clarity to the development team who are used to
    // the semantics of this function -CLG
    makeItSo();
}

template<class OP>
CopierIterator
Copier<OP>::begin(MotionType a_type) const
{
    switch (a_type)
    {
        case FROM:
            return CopierIterator(&m_fromMotionPlan);
        case TO:
            return CopierIterator(&m_toMotionPlan);
        case LOCAL:
            return CopierIterator(&m_localMotionPlan);
        default:
            PROTO_ASSERT(false, "Copier::begin | Error: Invalid MotionType");
    }
}
//========================================================================
// COPIER PRIVATE FUNCTIONS

template<class OP>
void
Copier<OP>::makeItSo()
{
    makeItSoBegin();
    makeItSoLocal();
    makeItSoEnd();
}

template<class OP>
void
Copier<OP>::makeItSoBegin()
{
#ifdef PR_MPI
    allocateBuffers();
    writeToSendBuffers();
    {
        PR_TIME("post_messages");
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
    }
#endif
}

template<class OP>
void
Copier<OP>::makeItSoLocal()
{
    PR_TIME("Copier::makeItSoLocal");
    for (auto iter = begin(LOCAL); iter.ok(); ++iter)
    {
        const auto& item = *iter;
        m_op.localCopy(item.fromRegion, item.fromIndex, item.toRegion, item.toIndex);
    }
}

template<class OP>
void
Copier<OP>::makeItSoEnd()
{
#ifdef PR_MPI
   completePendingSends();
   readFromRecvBuffers(); 
#endif
}

template<class OP>
void
Copier<OP>::completePendingSends() const
{
#ifdef PR_MPI
    PR_TIME("Copier::completePendingSends");
    if (m_numSends > 0)
    {
        PR_TIME("MPI_Waitall");
        m_sendStatus.resize(m_numSends);
        int result = MPI_Waitall(m_numSends, &m_sendRequests[0], &m_sendStats[0]);
        if (result != MPI_SUCCESS)
        {
            // TODO: Figure out what to do here. 
        }
    }
    m_numSends = 0;
#endif
}

template<class OP>
void
Copier<OP>::allocateBuffers()
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
        BufferEntry b;
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
        BufferEntry b;
        b.item = &item;
        b.size = m_op.linearSize(item.fromRegion, item.fromIndex);
        recvBufferSize += b.size;
        b.procID = item.procID;
        m_toMe.push_back(b);
    }
    std::sort(m_toMe.begin(); m_toMe.end());

    // allocate send and recv buffers if needed
    if (sendBufferSize > m_sendCapacity)
        //free(m_sendBuffer);
        //m_sendBuffer = malloc(sendBufferSize);
        m_sendBuffer = std::make_unique<char[]>(sendBufferSize);
        if (m_sendBuffer == NULL)
        {
            MayDay<void>::Error("Copier::allocateBuffers | Error: Out of memory.");
        }
        m_sendCapacity = sendBufferSize;
    }
    if (recvBufferSize > m_recvCapacity)
        //free(m_recvBuffer);
        //m_recvBuffer = malloc(recvBufferSize);
        m_recvBuffer = std::make_unique<char[]>(recvBufferSize);
        if (m_recvBuffer == NULL)
        {
            MayDay<void>::Error("Copier::allocateBuffers | Error: Out of memory.");
        }
        m_recvCapacity = recvBufferSize;
    }
#endif
}

template<class OP>
void
Copier<OP>::writeToSendBuffers()
{
#ifdef PR_MPI
    PR_TIME("Copier::writeToSendBuffers");
    for (unsigned int ii = 0; ii < m_fromMe.size(); ii++)
    {
        auto item = m_fromMe[ii];
        m_op.linearOut(item.bufferPointer,
            item.item->fromRegion,
            item.item->fromIndex);
    }
#endif
}

template<class OP>
void
Copier<OP>::postSends() const
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
        char* buffer = (char*)entry.bufferPointer;
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

template<class OP>
void
Copier<OP>::readFromRecvBuffers()
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
            m_op.linearIn(entry.bufferPointer,
                entry.item->toRegion,
                entry.item->toIndex);
        }
    }
    m_numRecvs = 0;
#endif
}

template<class OP>
void
Copier<OP>::postRecvs() const
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
        char* buffer = (char*)entry.bufferPointer;
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
            //maxSize = PR_MAX_MPI_MESSAGE_SIZE;
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
        //maxSize = std::max<long long>(bufferSize, maxSize);
        // skip all the size==0 messages from coalescence
        while(next < m_toMe.size() && m_toMe[next].size == 0) { ++next; }
    }
    for (auto iter = extraRequests.begin(); iter != extraRequests.end(); ++iter)
    {
        m_recvRequests.push_back(*iter);
    }
    m_numRecvs = m_recvRequests.size();
#endif
}

//========================================================================
// MOTIONITEM, BUFFERENTRY, and COPYITERATOR FUNCTIONS

MotionItem::MotionItem(const DataIndex& a_fromIndex,
        const DataIndex& a_toIndex,
        const Box& a_fromRegion,
        const Box& a_toRegion)
{
    fromIndex  = a_fromIndex;
    toIndex    = a_toIndex;
    fromRegion = a_fromRegion;
    toRegion   = a_toRegion;
    procID     = -1;
}

bool
MotionItem::operator==(const MotionItem& a_rhs) const
{
    return (fromIndex  == a_rhs.fromIndex) &&
           (toIndex    == a_rhs.toIndex) &&
           (fromRegion == a_rhs.fromRegion) &&
           (toRegion   == a_rhs.toRegion) &&
           (procID     == a_rhs.procID);
}

bool
MotionItem::operator<(const MotionItem& a_rhs) const
{
    if (procID == a_rhs.procID)
    {
        auto left = toRegion;
        auto right = a_rhs.toRegion;
        if (left.low() < right.low())
        {
            return (left.high() < right.high());
        } else {
            return (left < right);
        }
    } else {
        return (procID < a_rhs.procID);
    }
}


bool
BufferEntry::operator<(const BufferEntry& a_rhs) const
{
    return ((*item) < *(a_rhs.item));
}


CopierIterator::CopierIterator(const std::vector<MotionItem*>* a_motionPlan)
{
    m_current = 0;
    m_motionPlan = a_motionPlan;
}

const MotionItem&
CopierIterator::operator*() const
{
    PROTO_ASSERT(m_current < size(),
        "CopierIterator::operator() | Error: Index out of bounds: %i >= %lu.",
        m_current, size());
    return *(m_motionPlan->operator[](m_current));
}

const MotionItem&
CopierIterator::operator()() const
{
    return *(*this);
}

const MotionItem&
CopierIterator::operator[](size_t a_index) const
{
    PROTO_ASSERT(a_index < size(),
        "CopierIterator::operator() | Error: Index out of bounds: %lu >= %lu.",
        a_index, size());
    return *(m_motionPlan->operator[](a_index));
}

void
CopierIterator::operator++()
{
    ++m_current;
}

bool
CopierIterator::ok()
{
    return (m_current < size());
}

void
CopierIterator::reset()
{
    m_current = 0;
}

size_t
CopierIterator::size() const
{
    return m_motionPlan->size();
}


