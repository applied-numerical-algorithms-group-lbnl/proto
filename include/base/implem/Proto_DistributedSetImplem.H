template <typename DATA>
DistributedSet<DATA>::DistributedSet()
{
#ifdef PR_MPI
    m_sendBufferSizeInBytes = 0;
    m_recvBufferSizeInBytes = 0;
    m_sendBuffer = nullptr;
    m_recvBuffer = nullptr;
    protoMalloc(MEMTYPE_DEFAULT, m_localBufferSizesInBytes, numProc() * sizeof(int));
    protoMalloc(MEMTYPE_DEFAULT, m_bufferOffsetsInBytes, numProc() * sizeof(int));
#endif
}

template <typename DATA>
DistributedSet<DATA>::~DistributedSet()
{
    proto_free(m_localBufferSizesInBytes);
    proto_free(m_bufferOffsetsInBytes);
    if (m_sendBufferSizeInBytes > 0)
    {
        proto_free(m_sendBuffer);
    }
    if (m_recvBufferSizeInBytes > 0)
    {
        proto_free(m_recvBuffer);
    }
}

template <typename DATA>
void DistributedSet<DATA>::insert(const DATA &item)
{
#ifdef PR_MPI
    m_localBuffer.push_back(item);
#else
    m_globalSet.insert(item);
#endif
    return;
}

template <typename DATA>
void DistributedSet<DATA>::exchange()
{
#ifdef PR_MPI
    if (numProc() == 1)
    {
        for (auto item : m_localBuffer)
        {
            m_globalSet.insert(item);
        }
        return;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // GET BUFFER SIZES
    int localBufferSizeInBytes = m_localBuffer.size() * sizeof(DATA);
    MPI_Allgather(&localBufferSizeInBytes, 1, MPI_INT, m_localBufferSizesInBytes, 1, MPI_INT, MPI_COMM_WORLD);

    // ALLOCATE SEND BUFFER
    if (localBufferSizeInBytes > m_sendBufferSizeInBytes)
    {
        if (m_sendBufferSizeInBytes > 0)
        {
            protoFree(MEMTYPE_DEFAULT, m_sendBuffer);
        }
        m_sendBufferSizeInBytes = localBufferSizeInBytes;
        protoMalloc(MEMTYPE_DEFAULT, m_sendBuffer, m_sendBufferSizeInBytes);
    }

    // WRITE TO SEND BUFFER
    int writeIndex = 0;
    for (auto item : m_localBuffer)
    {
        proto_memcpy(&item, m_sendBuffer + writeIndex, sizeof(DATA));
        writeIndex += sizeof(DATA);
    }

    // ALLOCATE GLOBAL BUFFER AND OFFSETS
    int globalBufferSizeInBytes = 0;
    for (int ii = 0; ii < numProc(); ii++)
    {
        m_bufferOffsetsInBytes[ii] = globalBufferSizeInBytes;
        globalBufferSizeInBytes += m_localBufferSizesInBytes[ii];
    }
    if (globalBufferSizeInBytes > m_recvBufferSizeInBytes)
    {
        if (m_recvBufferSizeInBytes > 0)
        {
            protoFree(MEMTYPE_DEFAULT, m_recvBuffer);
        }
        m_recvBufferSizeInBytes = globalBufferSizeInBytes;
        protoMalloc(MEMTYPE_DEFAULT, m_recvBuffer, m_recvBufferSizeInBytes);
    }

    // EXCHANGE DATA
    MPI_Allgatherv(m_sendBuffer, localBufferSizeInBytes, MPI_BYTE,
                   m_recvBuffer, m_localBufferSizesInBytes, m_bufferOffsetsInBytes, MPI_BYTE, MPI_COMM_WORLD);

    // READ DATA FROM BUFFER
    int readIndex = 0;
    for (int readIndex = 0; readIndex < globalBufferSizeInBytes; readIndex += sizeof(DATA))
    {
        DATA item;
        proto_memcpy(m_recvBuffer + readIndex, &item, sizeof(DATA));
        m_globalSet.insert(item);
    }
#endif
}

template <typename DATA>
void DistributedSet<DATA>::clearLocal()
{
#ifdef PR_MPI
    m_localBuffer.clear();
#endif
}

template <typename DATA>
void DistributedSet<DATA>::clearGlobal()
{
    m_globalSet.clear();
}