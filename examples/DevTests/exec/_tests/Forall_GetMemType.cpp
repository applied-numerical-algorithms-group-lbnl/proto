#include "Proto.H"
#include "InputParser.H"

using namespace Proto;

template<typename T>
MemType get(T& a_dummy)
{
    return getMemType<T>::type_eval();
}

int main(int argc, char** argv)
{
    #ifdef PR_MPI
    MPI_Init(&argc, &argv);
    #endif

    BoxData<double, 1, HOST> hostData;
    BoxData<double, 1, DEVICE> deviceData;
    Var<double, 1, HOST> hostVar;
    Var<double, 1, DEVICE> deviceVar;

    MemType mem_hostData    = getMemType<BoxData<double, 1, HOST, 1, 1>>::type_eval();
    MemType mem_deviceData  = getMemType<BoxData<double, 1, DEVICE, 1, 1>>::type_eval();
    MemType mem_hostVar     = getMemType<decltype(hostVar)>::type_eval();
    MemType mem_deviceVar   = getMemType<decltype(deviceVar)>::type_eval();
    MemType mem_double      = getMemType<double>::type_eval();
    MemType mem_int         = getMemType<int>::type_eval();
    MemType mem_Point       = getMemType<Point>::type_eval();
    MemType mem_Box         = getMemType<Box>::type_eval();
    MemType mem_array       = getMemType<std::array<double, DIM>>::type_eval();
    
    bool pass = true;
    pass &= (mem_hostData   == HOST);
    pass &= (mem_deviceData == DEVICE);
    pass &= (mem_hostVar    == HOST);
    pass &= (mem_deviceVar  == DEVICE);
    pass &= (mem_double     == BOTH );
    pass &= (mem_int        == BOTH );
    pass &= (mem_Point      == BOTH );
    pass &= (mem_Box        == BOTH );
    pass &= (mem_array      == BOTH );

    if (procID() == 0)
    {
        if (pass){std::cout << "ALL TESTS PASSED" << std::endl;}
        else 
        {
            std::cout << "TESTS FAILED" << std::endl;
            std::cout << "Host Data: "   << mem_hostData << std::endl;
            std::cout << "Device Data: " << mem_deviceData << std::endl;
            std::cout << "Host Var: "    << mem_hostVar << std::endl;
            std::cout << "Device Var: "  << mem_deviceVar << std::endl;
            std::cout << "double: "      << mem_double << std::endl;
            std::cout << "int: "         << mem_int << std::endl;
            std::cout << "Point: "       << mem_Point << std::endl;
            std::cout << "Box: "         << mem_Box << std::endl;
            std::cout << "array: "       << mem_array << std::endl;
        }
    }

    #ifdef PR_MPI
    MPI_Finalize();
#endif

    return 0;
}

