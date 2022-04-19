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
    
#if 0
    Reduction<int, Sum, DEVICE> rxn;
    int sum = 0;
    rxn.reduce(data.data(), data.size());
    sum = rxn.update(sum, rxn.fetch());
    rxn.reset();
    std::cout << "Sum: " << sum << " | Should be: " << domainBox.size() << std::endl;
    rxn.reduce(data.data(), data.size());
    sum = rxn.update(sum, rxn.fetch());
    rxn.reset();
    std::cout << "Sum: " << sum << " | Should be: " << domainBox.size()*2 << std::endl;
    rxn.reduce(data.data(), data.size());
    sum = rxn.update(sum, rxn.fetch());
    rxn.reset();
    std::cout << "Sum: " << sum << " | Should be: " << domainBox.size()*3 << std::endl;
    rxn.reduce(data.data(), data.size());
    sum = rxn.update(sum, rxn.fetch());
    rxn.reset();
    std::cout << "Sum: " << sum << " | Should be: " << domainBox.size()*4 << std::endl;
#else
    Reduction<int, Sum, DEVICE> rxn;
    rxn.reduce(data.data(), data.size());
    std::cout << "Sum: " << rxn.fetch() << " | Should be: " << domainBox.size() << std::endl;
    rxn.reduce(data.data(), data.size());
    std::cout << "Sum: " << rxn.fetch() << " | Should be: " << domainBox.size()*2 << std::endl;
    rxn.reduce(data.data(), data.size());
    std::cout << "Sum: " << rxn.fetch() << " | Should be: " << domainBox.size()*3 << std::endl;
    rxn.reduce(data.data(), data.size());
    std::cout << "Sum: " << rxn.fetch() << " | Should be: " << domainBox.size()*4 << std::endl;
#endif
    #ifdef PR_MPI
    MPI_Finalize();
    #endif
    return 0;
}

