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
        std::cout << "usage: ./testLevelData domainSize boxSize" << std::endl;
        if (argc > 1)
        {
            domainSize = atoi(argv[1]);
        }
        if (argc > 2)
        {
            boxSize = atoi(argv[2]);
        }
        std::cout << "domainSize = " << domainSize << " | boxSize = " << boxSize << std::endl;
    }

#ifdef PR_MPI
    MPI_Bcast(&domainSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&boxSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

    HDF5Handler h5;

    Box domain = Box::Cube(domainSize);
    Point boxSizeVect = Point::Ones(boxSize);

    Point ghost = Point::Ones();
    std::array<bool, DIM> periodicity;
    for (int ii = 0; ii < DIM; ii++) { periodicity[ii] = true; }
    
    ProblemDomain problemDomain(domain, periodicity);
    DisjointBoxLayout layout(problemDomain, boxSizeVect);
    
    LevelBoxData<double> src(layout, ghost);
    LevelBoxData<double> dst(layout, ghost);
    for (auto iter = src.begin(); iter.ok(); ++iter)
    {
        auto& patch = src[*iter];
        patch.setVal(*iter, layout[*iter]);
    }

    h5.writeLevel(src, "SRC_0");
    h5.writeLevel(dst, "DST_0");
    src.exchange();
    h5.writeLevel(src, "SRC_1");
    src.copyTo(dst);
    h5.writeLevel(dst, "DST_1");

    bool success = true;
    Reduction<double> rxn;
    for (auto iter = layout.begin(); iter.ok(); ++iter)
    {
        auto& patch = src[*iter];
        patch.absMax(rxn);
        //double max = rxn.fetch();
        double max = patch.absMax();
        std::cout << "Proc: " << procID() << " | patch: " << *iter << " | abs-max: " << max << std::endl;
    }
#ifdef PR_MPI
    MPI_Finalize();
#endif
}
