#include "Proto.H"
#include "InputParser.H"

using namespace Proto;

int main(int argc, char** argv)
{
    #ifdef PR_MPI
    MPI_Init(&argc, &argv);
    #endif

    // SETUP
    HDF5Handler h5;

    int domainSize = 8;
    InputArgs args;
    args.add("domainSize", domainSize);
    args.parse(argc, argv);
    args.print();

    Box domainBox = Box::Cube(domainSize);
    BoxData<int, 1, DEVICE> data(domainBox);
    data.setVal(1);

    Reduction<int, Sum, DEVICE> rxn;

    rxn.reduce(data.data(), data.linearSize());
    
    std::cout << "Sum: " << rxn.fetch() << " | Should be: " << domainBox.size() << std::endl;

    #ifdef PR_MPI
    MPI_Finalize();
    #endif
    return 0;
}

