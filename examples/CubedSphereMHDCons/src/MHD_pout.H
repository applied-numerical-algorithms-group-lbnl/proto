#ifdef PR_MPI
#include <mpi.h>
#endif

// ======================================================================
//  Per-rank pout() implementation (global namespace)
// ======================================================================

// Null stream buffer that discards everything
class HelioCubedNullBuffer : public std::streambuf
{
protected:
    int overflow(int c) override
    {
        // Pretend we successfully wrote the character
        return traits_type::not_eof(c);
    }
};

// Base pout(): one file per rank, or std::cout if no MPI / not initialized
inline std::ostream& pout()
{
    static thread_local std::ostream* os = nullptr;

    if (os != nullptr)
    {
        return *os;
    }

#ifdef PR_MPI
    int mpi_initialized = 0;
    MPI_Initialized(&mpi_initialized);

    if (mpi_initialized)
    {
        int rank = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        std::ostringstream fname;
        fname << "pout." << rank;

        // New file each run: truncate instead of append
        auto ofs = new std::ofstream(fname.str(),
                                     std::ios::out | std::ios::trunc);

        if (!ofs->is_open())
        {
            // On failure, fall back to stdout
            std::cerr << "WARNING: could not open " << fname.str()
                      << ", falling back to std::cout\n";
            delete ofs;
            os = &std::cout;
        }
        else
        {
            ofs->setf(std::ios::unitbuf); // flush on each insertion
            os = ofs;
        }

        return *os;
    }
#endif
    // Non-MPI build, or MPI not initialized: use stdout
    os = &std::cout;
    return *os;
}

// pout(rootRank): only rootRank writes; others discard output
inline std::ostream& pout(int rootRank)
{
#ifdef PR_MPI
    int mpi_initialized = 0;
    MPI_Initialized(&mpi_initialized);

    if (mpi_initialized)
    {
        int rank = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        if (rank != rootRank)
        {
            // Non-root ranks write to a null stream
            static thread_local HelioCubedNullBuffer nullBuf;
            static thread_local std::ostream nullStream(&nullBuf);
            return nullStream;
        }

        // This is the selected rank: write to its regular pout.<rank> file
        return pout();
    }
#endif
    // No MPI (or not initialized): treat as single process = rank 0
    if (rootRank == 0)
    {
        return pout(); // whatever pout() is using (usually std::cout)
    }
    else
    {
        static thread_local HelioCubedNullBuffer nullBuf;
        static thread_local std::ostream nullStream(&nullBuf);
        return nullStream;
    }
}