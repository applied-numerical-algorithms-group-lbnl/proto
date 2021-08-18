#include "Proto.H"

using namespace Proto;

int main(int argc, char** argv)
{
#ifdef PR_MPI
    MPI_Init(&argc, &argv);
#endif
    int domainSize = 16;
    Box domain = Box::Cube(domainSize);
    Point boxSize = Point::Ones(4);
    Point ghost = Point::Ones();
    std::array<bool, DIM> periodicity;
    for (int ii = 0; ii < DIM; ii++) { periodicity[ii] = true; }
    
    ProblemDomain problemDomain(domain, periodicity);
    DisjointBoxLayout layout(problemDomain, boxSize);

    for (auto iter = layout.begin(); iter.ok(); ++iter)
    {
        if (procID() == 0)
        {
            auto b0 = layout[*iter];
            auto b =  layout.domain() & b0.grow(ghost);
            std::cout << "Source Box = " << b0 << std::endl;
            std::cout << "Neighbors: " << std::endl;
            NeighborIterator niter(layout, b);
            for (niter.begin(); *niter != niter.end(); ++niter)
            {
                std::cout << "src box: " << niter.srcBox() << " | dst box: " << niter.destBox() << " | shift: " << niter.shift() << " | proc: " << niter.procID() << std::endl;
            }
        }
    }

    //LevelBoxData<double> src(layout, ghost);
    //LevelBoxData<double> dst(layout, ghost);

    for (auto iter = layout.begin(); iter.ok(); ++iter)
    {
        //src[*iter].setVal((*iter)+1, layout[*iter]);
        //dst[*iter].setVal(0);
    }

    //HDF5Handler h5;

    //h5.writeLevel(src, "SRC_0.hdf5");
    //h5.writeLevel(dst, "DST_0.hdf5");

    //src.exchange();
    //src.copyTo(dst);
    
    //h5.writeLevel(src, "SRC_1.hdf5");
    //h5.writeLevel(dst, "DST_1.hdf5");

#ifdef PR_MPI
    MPI_Finalize();
#endif
}
