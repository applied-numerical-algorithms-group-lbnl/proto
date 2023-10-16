#include "Proto.H"

int main(int argc, char** argv)
{
    #ifdef PR_MPI
    MPI_Init(&argc, &argv);
    #endif
    
    pr_out() << "Hello from rank " << procID() << std::endl;

#ifdef PR_MPI
    MPI_Finalize();
    #endif
    return 0;
}
