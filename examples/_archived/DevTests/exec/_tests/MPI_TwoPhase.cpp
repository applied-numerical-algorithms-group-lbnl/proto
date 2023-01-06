#include "Proto.H"
#include "InputParser.H"

using namespace Proto;

int main(int argc, char** argv)
{


    // m_bufferSizes[numProc]
    // m_localBuffer;
    // m_globalBuffer;
    //

    #ifdef PR_MPI
    MPI_Init(&argc, &argv);
    
    int proc = procID();
    int bufferSizes[numProc()];
    int bufferSize = proc;
    
    MPI_Allgather(&bufferSize, 1, MPI_INT, bufferSizes, 1, MPI_INT, MPI_COMM_WORLD);

    int localBuffer[bufferSize];
    for (int ii = 0; ii < bufferSize; ii++)
    {
        localBuffer[ii] = proc;
    }

    int totalBufferSize = 0;
    int offsets[numProc()];
    for (int ii = 0; ii < numProc(); ii++)
    {
        offsets[ii] = totalBufferSize;
        totalBufferSize += bufferSizes[ii];
    }

    int totalBuffer[totalBufferSize];

    MPI_Allgatherv(localBuffer, bufferSize, MPI_INT, totalBuffer, bufferSizes, offsets, MPI_INT, MPI_COMM_WORLD);

    pout() << "Final Buffer" << std::endl;
    for (int ii = 0; ii < totalBufferSize; ii++)
    {
        pout() << totalBuffer[ii] << ", ";
    }
    pout() << std::endl;


    MPI_Finalize();
    #endif
}
