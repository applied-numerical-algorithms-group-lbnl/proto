#include "Proto.H"
#include "InputParser.H"
#include "TestFunc.H"

using namespace Proto;

int main(int argc, char** argv)
{
    #ifdef PR_MPI
    MPI_Init(&argc, &argv);
    #endif

    // SETUP
    HDF5Handler h5;

    int domainSize = 64;
    double physDomainSize = 1.0;
    int testNum = 0;
    int numIter = 3;
    InputArgs args;
    args.add("domainSize",      domainSize);
    args.add("physDomainSize",  physDomainSize);
    args.add("testNum",         testNum);
    args.add("numIter",         numIter);
    args.parse(argc, argv);
    args.print();

    if (testNum == 0 && procID() == 0)
    {
        std::cout << "TEST 0: MEMTYPE = HOST" << std::endl;
    }
    else if (testNum == 1 && procID() == 0)
    {
        std::cout << "TEST 1: MEMTYPE = DEVICE" << std::endl;
    }

    Stencil<double> L = Stencil<double>::Laplacian();
    double k = 2.0*M_PI;

    double err_replace[numIter];
    double err_add[numIter];
    double err_assign[numIter];
    for (int nn = 0; nn < numIter; nn++)
    {
        double dx = physDomainSize / domainSize;
        Box domainBox = Box::Cube(domainSize);
        if (testNum == 0)
        {
            BoxData<double, 1, HOST> src(domainBox.grow(1));
            BoxData<double, 1, HOST> dst_replace(domainBox.grow(1));
            BoxData<double, 1, HOST> dst_add(domainBox.grow(1));
            BoxData<double, 1, HOST> sln_replace(domainBox.grow(1));
            BoxData<double, 1, HOST> sln_add(domainBox.grow(1));
            BoxData<double, 1, HOST> sln_assign(domainBox);

            forallInPlace_p(sinProd_avg, src, dx);
            sln_assign.setVal(0);
            sln_assign += src;
            sln_add.setVal(7);
            sln_add += sln_assign;
            sln_replace.setVal(7);
            forallInPlace_p(sinProd_avg, domainBox, sln_replace, dx);

            dst_add.setVal(7);
            dst_replace.setVal(7);

            BoxData<double, 1, HOST> dst_assign = L(src, -1/(DIM*dx*dx*k*k));
            dst_replace |= L(src, -1/(DIM*dx*dx*k*k)); 
            dst_add += L(src, -1/(DIM*dx*dx*k*k)); 

            dst_assign -= sln_assign;
            dst_replace -= sln_replace;
            dst_add -= sln_add;
            
            err_assign[nn] = dst_assign.absMax();
            err_add[nn] = dst_add.absMax();
            err_replace[nn] = dst_replace.absMax();

            if (procID() == 0)
            {
                std::cout << "assign test error: " << err_assign[nn] << std::endl;
                std::cout << "replace test error: " << err_replace[nn] << std::endl;
                std::cout << "add test error: " << err_add[nn] << std::endl;
            }

            bool PASS = true;
            PASS &= (sln_assign.box() == dst_assign.box());
            if (!PASS && procID() == 0)
            {
                std::cout << "FAILURE | assignment box size test failed | correct box: " << domainBox << " | computed box: " << dst_assign.box() << std::endl;
            }
        }
        else if (testNum == 1)
        {
#ifdef PROTO_ACCEL
            BoxData<double, 1, DEVICE> src(domainBox.grow(1));
            BoxData<double, 1, DEVICE> dst_replace(domainBox.grow(1));
            BoxData<double, 1, DEVICE> dst_add(domainBox.grow(1));
            BoxData<double, 1, DEVICE> sln_replace(domainBox.grow(1));
            BoxData<double, 1, DEVICE> sln_add(domainBox.grow(1));
            BoxData<double, 1, DEVICE> sln_assign(domainBox);

            forallInPlace_p(sinProd_avg, src, dx);
            sln_assign.setVal(0);
            sln_assign += src;
            sln_add.setVal(7);
            sln_add += sln_assign;
            sln_replace.setVal(7);
            forallInPlace_p(sinProd_avg, domainBox, sln_replace, dx);

            dst_add.setVal(7);
            dst_replace.setVal(7);

            BoxData<double, 1, DEVICE> dst_assign = L(src, -1/(DIM*dx*dx*k*k));
            dst_replace |= L(src, -1/(DIM*dx*dx*k*k)); 
            dst_add += L(src, -1/(DIM*dx*dx*k*k)); 

            dst_assign -= sln_assign;
            dst_replace -= sln_replace;
            dst_add -= sln_add;
            
            err_assign[nn] = dst_assign.absMax();
            err_add[nn] = dst_add.absMax();
            err_replace[nn] = dst_replace.absMax();

            if (procID() == 0)
            {
                std::cout << "assign test error: " << err_assign[nn] << std::endl;
                std::cout << "replace test error: " << err_replace[nn] << std::endl;
                std::cout << "add test error: " << err_add[nn] << std::endl;
            }

            bool PASS = true;
            PASS &= (sln_assign.box() == dst_assign.box());
            if (!PASS && procID() == 0)
            {
                std::cout << "FAILURE | assignment box size test failed | correct box: " << domainBox << " | computed box: " << dst_assign.box() << std::endl;
            }
#else
            std::cout << "Test not run: CUDA is not enabled." << std::endl;
#endif
        }
        domainSize *= 2;
    }
    
    for (int ii = 1; ii < numIter; ii++)
    {
        if (procID() == 0)
        {
            std::cout << "Convergence Rate (Assign): " << log(err_assign[ii-1] / err_assign[ii]) / log(2.0) << std::endl;
            std::cout << "Convergence Rate (Replace): " << log(err_replace[ii-1] / err_replace[ii]) / log(2.0) << std::endl;
            std::cout << "Convergence Rate (Add): " << log(err_add[ii-1] / err_add[ii]) / log(2.0) << std::endl;
        }
    }
    #ifdef PR_MPI
    MPI_Finalize();
    #endif
    return 0;
}

