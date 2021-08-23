#include "Proto.H"

using namespace Proto;

int main(int argc, char** argv)
{
#ifdef PR_MPI
    MPI_Init(&argc, &argv);
#endif
    int domainSize;
    int boxSize;
    if (procID() == 0)
    {
        if (argc > 1)
        {
            domainSize = atoi(argv[1]);
        }
        if (argc > 2)
        {
            boxSize = atoi(argv[2]);
        }
    }

#ifdef PR_MPI
    MPI_Bcast(&domainSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&boxSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

    Box domain = Box::Cube(domainSize);
    Point boxSizeVect = Point::Ones(boxSize);

    Point ghost = Point::Ones();
    std::array<bool, DIM> periodicity;
    for (int ii = 0; ii < DIM; ii++) { periodicity[ii] = true; }
    
    ProblemDomain problemDomain(domain, periodicity);
    DisjointBoxLayout layout(problemDomain, boxSizeVect);



#ifdef PR_MPI
    MPI_Finalize();
#endif
}
