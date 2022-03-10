
FinitePointSet::FinitePointSet(Box a_domain, std::array<bool, DIM> a_periodic)
{
    m_domain = ProblemDomain(a_domain, a_periodic);
#ifdef PR_MPI
    m_sendBufferSize = 0;
    m_recvBufferSize = 0;
    m_sendBuffer = nullptr;
    m_recvBuffer = nullptr;
    
    m_globalBufferSize = 0;
    m_localBufferSize = 0;
    m_bufferSizes = (int*)malloc(numProc()*sizeof(int));
    m_bufferOffsets = (int*)malloc(numProc()*sizeof(int));
    m_localBuffer = nullptr;
    m_globalBuffer = nullptr;
#endif
}

FinitePointSet::FinitePointSet(Box a_domain, bool a_periodic)
{
    m_domain = ProblemDomain(a_domain, a_periodic);
#ifdef PR_MPI
    m_sendBufferSize = 0;
    m_recvBufferSize = 0;
    m_sendBuffer = nullptr;
    m_recvBuffer = nullptr;
    
    m_globalBufferSize = 0;
    m_localBufferSize = 0;
    m_bufferSizes = (int*)malloc(numProc()*sizeof(int));
    m_bufferOffsets = (int*)malloc(numProc()*sizeof(int));
    m_localBuffer = nullptr;
    m_globalBuffer = nullptr;
#endif
}

FinitePointSet::FinitePointSet(ProblemDomain a_domain)
{
    m_domain = a_domain;
#ifdef PR_MPI
    m_sendBufferSize = 0;
    m_recvBufferSize = 0;
    m_sendBuffer = nullptr;
    m_recvBuffer = nullptr;
    
    m_globalBufferSize = 0;
    m_localBufferSize = 0;
    m_bufferSizes = (int*)malloc(numProc()*sizeof(int));
    m_bufferOffsets = (int*)malloc(numProc()*sizeof(int));
    m_localBuffer = nullptr;
    m_globalBuffer = nullptr;
#endif
}

FinitePointSet::~FinitePointSet()
{
#ifdef PR_MPI
    if (m_sendBuffer != nullptr)
    {
        free(m_sendBuffer);
    }
    if (m_recvBuffer != nullptr)
    {
        free(m_recvBuffer);
    }
    free(m_bufferSizes);
    free(m_bufferOffsets);
    if (m_localBufferSize > 0)
    {
        free(m_localBuffer);
    }
    if (m_globalBufferSize > 0)
    {
        free(m_globalBuffer);
    }
#endif
}

bool FinitePointSet::add(Point a_point)
{
    if (m_points.count(a_point) > 0) {return false; }
    Point p = m_domain.image(a_point); 
    if (!m_domain.box().contains(p))
    {
        return false;
    }
#ifdef PR_MPI
    m_pointBuffer.push_back(a_point);
#else
    m_points.insert(a_point);
#endif
    return true;
}

void FinitePointSet::exchange()
{
#ifdef PR_MPI
    MPI_Barrier(MPI_COMM_WORLD);
    m_points.clear();

    // GET BUFFER SIZES
    int localBufferSize = DIM*m_pointBuffer.size();
    MPI_Allgather(&localBufferSize, 1, MPI_INT, m_bufferSizes, 1, MPI_INT, MPI_COMM_WORLD);  

    // ALLOCATE LOCAL BUFFER
    if (localBufferSize > m_localBufferSize)
    {
        if (m_localBufferSize > 0)
        {
            free(m_sendBuffer);
        }
        m_localBufferSize = localBufferSize;
        m_localBuffer = (int*)malloc(m_localBufferSize*sizeof(int));
    }

    // WRITE TO LOCAL BUFFER
    int writeIndex = 0;
    for (auto iter = m_pointBuffer.begin(); iter != m_pointBuffer.end(); ++iter)
    {
        char* buf = ((char*)m_localBuffer) + writeIndex;
        (*iter).linearOut(buf);
        writeIndex += DIM*sizeof(int);
    }

    // ALLOCATE GLOBAL BUFFER AND OFFSETS
    int globalBufferSize = 0;
    for (int ii = 0; ii < numProc(); ii++)
    {
        m_bufferOffsets[ii] = globalBufferSize;
        globalBufferSize += m_bufferSizes[ii];
    }
    if (globalBufferSize > m_globalBufferSize)
    {
        if (m_globalBufferSize > 0)
        {
            free(m_globalBuffer);
        }
        m_globalBufferSize = globalBufferSize;
        m_globalBuffer = (int*)malloc(m_globalBufferSize*sizeof(int));
    }

    // EXCHANGE DATA
    MPI_Allgatherv(m_localBuffer, localBufferSize, MPI_INT, 
            m_globalBuffer, m_bufferSizes, m_bufferOffsets, MPI_INT, MPI_COMM_WORLD);

    // READ DATA FROM BUFFER
    int readIndex = 0;
    for (int ii = 0; ii < globalBufferSize / DIM; ii++)
    {
        char* buf = ((char*)m_globalBuffer) + readIndex;
        Point p;
        p.linearIn(buf);
        m_points.insert(p);
        readIndex += DIM*sizeof(int);
    }
#endif
}

void FinitePointSet::exchange_old()
{
#ifdef PR_MPI
    MPI_Barrier(MPI_COMM_WORLD);
    // LOCAL INSERTIONS
    pout() << "============================================================" << std::endl;
    pout() << "FinitePointSet::exchange begin" << std::endl;
    pout() << "Sending " << m_pointBuffer.size() << " local points:" << std::endl;
    for (int ii = 0; ii < m_pointBuffer.size(); ii++)
    {
        pout() << "\t" << m_pointBuffer[ii] << std::endl;
        m_points.insert(m_pointBuffer[ii]);
    }
    
    // ALLOC SEND BUFFER (IF NEEDED)
    int sendBufferSize = DIM*m_pointBuffer.size();
    if (sendBufferSize > m_sendBufferSize)
    {
        if (m_sendBufferSize > 0)
        {
            free(m_sendBuffer);
        }
        m_sendBufferSize = sendBufferSize;
        m_sendBuffer = (int*)malloc(m_sendBufferSize*sizeof(int));
    }
    // FILL SEND BUFFER
    int* pos = m_sendBuffer;
    for (int ii = 0; ii < m_pointBuffer.size(); ii++)
    {
        Point& p = m_pointBuffer[ii];
        for (int d = 0; d < DIM; d++)
        {
            pos[d] = p[d];
        }
        pos += DIM;
    }
    
    // LAUNCH ASYNCHRONOUS SENDS
    std::vector<MPI_Request> sendRequests(Proto::numProc());
    for (int ii = 0; ii < Proto::numProc(); ii++)
    {
        if (ii == Proto::procID()) {continue;}
        pout() << "Sending data | count: " << sendBufferSize << std::endl;
        MPI_Isend((void*)m_sendBuffer, sendBufferSize, MPI_INT, 
                ii, ii, MPI_COMM_WORLD, &sendRequests[ii]);
    }

    // PROBE FOR RECVS
    std::vector<int> recvStatus(Proto::numProc(), false);
    recvStatus[Proto::procID()] = true;
    bool wait = true;
    MPI_Barrier(MPI_COMM_WORLD);
    pout() << "Waiting for data..." << std::endl;
    while (wait)
    {
        for (int ii = 0; ii < Proto::numProc(); ii++)
        {
            MPI_Status status;
            if (ii == Proto::procID()) {continue;}
            MPI_Iprobe(ii, procID(), MPI_COMM_WORLD, &recvStatus[ii], &status);
            if (recvStatus[ii])
            {
                int recvBufferSize = 0;
                MPI_Get_count(&status, MPI_INT, &recvBufferSize);
                pout() << "Found data | proc: " << ii << " | count: " << recvBufferSize << std::endl;
                if (recvBufferSize > m_recvBufferSize)
                {
                    if (m_recvBufferSize > 0)
                    {
                        free(m_recvBuffer);
                    }

                    m_recvBufferSize = recvBufferSize;
                    m_recvBuffer = (int*)malloc(recvBufferSize*sizeof(int));
                }
                if (recvBufferSize > 0)
                {
                    pout() << "\tReading Data" << std::endl;
                    MPI_Recv(m_recvBuffer, m_recvBufferSize, MPI_INT, ii, procID(),
                            MPI_COMM_WORLD, &status);
                    int numRecvPoints = recvBufferSize / DIM;
                    int* recvPos = m_recvBuffer;
                    for (int pi = 0; pi < numRecvPoints; pi++)
                    {
                        Point p;
                        for (int ii = 0; ii < DIM; ii++)
                        {
                            p[ii] = *(recvPos + ii);
                        }
                        m_points.insert(p);
                        recvPos += DIM; 
                    }
                }
            } 
        }
        pout() << "Are all procs done?";
        wait = false;
        for (int ii = 0; ii < Proto::numProc(); ii++)
        {
            if (ii == procID()) { continue; }
            if (recvStatus[ii] == 0)
            {
                pout() << " NO. Keep checking." << std::endl;
                wait = true;
                break;
            }
        }
        if (!wait) { pout() << " YES. Done checking." << std::endl;}
    }
    m_pointBuffer.clear();
#endif
}

std::vector<Point>
FinitePointSet::points()
{
    exchange();
    std::vector<Point> points;
    for (auto iter = m_points.begin(); iter != m_points.end(); ++iter)
    {
        points.push_back(*iter);
    }
    return points;
}

void FinitePointSet::clear()
{
    m_points.clear();
    m_pointBuffer.clear();
}

#ifdef PR_MPI
/*
int
FinitePointSet::recvBufferSize()
{
    // assume the largest message is twice the perfectly load balanced value
    return (m_domain.size() / Proto::numProc())*2*DIM;
}
*/
#endif

