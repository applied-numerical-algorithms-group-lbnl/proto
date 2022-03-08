
FinitePointSet::FinitePointSet(Box a_domain, std::array<bool, DIM> a_periodic)
{
    //m_domain = a_domain; 
    //m_periodic = a_periodic;
    m_domain = ProblemDomain(a_domain, a_periodic);
#ifdef PR_MPI
    m_sendBufferSize = 0;
    m_recvBufferSize = 0;
    m_sendBuffer = nullptr;
    m_recvBuffer = nullptr;
#endif
}

FinitePointSet::FinitePointSet(Box a_domain, bool a_periodic)
{
    //m_domain = a_domain;
    //for (int ii = 0; ii < DIM; ii++) { m_periodic[ii] = a_periodic; }
    m_domain = ProblemDomain(a_domain, a_periodic);
#ifdef PR_MPI
    m_sendBufferSize = 0;
    m_recvBufferSize = 0;
    m_sendBuffer = nullptr;
    m_recvBuffer = nullptr;
#endif
}

FinitePointSet::FinitePointSet(ProblemDomain a_domain)
{
    m_domain = a_domain;
    //for (int ii = 0; ii < DIM; ii++) { m_periodic[ii] = a_periodic; }
#ifdef PR_MPI
    m_sendBufferSize = 0;
    m_recvBufferSize = 0;
    m_sendBuffer = nullptr;
    m_recvBuffer = nullptr;
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
    // LOCAL INSERTIONS
    for (int ii = 0; ii < m_pointBuffer.size(); ii++)
    {
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
        MPI_Isend((void*)m_sendBuffer, sendBufferSize, MPI_INT, ii, ii, MPI_COMM_WORLD, &sendRequests[ii]);
    }

    // PROBE FOR RECVS
    std::vector<int> recvStatus(Proto::numProc(), false);
    recvStatus[Proto::procID()] = true;
    bool wait = true;
    MPI_Barrier(MPI_COMM_WORLD);
    while (wait)
    {
        for (int ii = 0; ii < Proto::numProc(); ii++)
        {
            MPI_Status status;
            if (ii == Proto::procID()) {continue;}
            MPI_Iprobe(ii, procID(), MPI_COMM_WORLD, &recvStatus[ii], &status);
            if (recvStatus[ii])
            {
                int recvBufferSize;
                MPI_Get_count(&status, MPI_INT, &recvBufferSize);
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
        wait = false;
        for (int ii = 0; ii < Proto::numProc(); ii++)
        {
            if (ii == procID()) { continue; }
            if (recvStatus[ii] == 0)
            {
                wait = true;
                break;
            }
        }
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

