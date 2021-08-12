#include "Proto.H"

using namespace Proto;

int main(int argc, char** argv)
{
#ifdef PR_MPI
    MPI_Init(&argc, &argv);
#endif
    int domainSize = 8;
    Box domain = Box::Cube(domainSize);
    Point boxSize = Point::Ones(4);
    Point ghost = Point::Ones();
    std::array<bool, DIM> periodicity;
    for (int ii = 0; ii < DIM; ii++) { periodicity[ii] = false; }
    
    ProblemDomain problemDomain(domain, periodicity);
    DisjointBoxLayout layout(problemDomain, boxSize);

    for (auto iter = layout.begin(); iter.ok(); ++iter)
    {
        std::cout << "process: " << procID() << layout[*iter] << std::endl;
    }

    LevelBoxData<double> src(layout, ghost);
    LevelBoxData<double> dst(layout, ghost);

    for (auto iter = layout.begin(); iter.ok(); ++iter)
    {
        src[*iter].setVal((*iter)+1, layout[*iter]);
        dst[*iter].setVal(0);
        src[*iter].printData();
    }
    src.exchange();
    for (auto iter = layout.begin(); iter.ok(); ++iter)
    {
        src[*iter].printData();
    }
    src.copyTo(dst);
    for (auto iter = layout.begin(); iter.ok(); ++iter)
    {
        dst[*iter].printData();
    }

#ifdef PR_MPI
    MPI_Finalize();
#endif
}
