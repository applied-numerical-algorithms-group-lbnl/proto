#include "ProtoAMR.H"
#include "InputParser.H"
#include "BoxOp_Laplace.H"
#include "TestFunc.H"

using namespace Proto;

int main(int argc, char** argv)
{
    #ifdef PR_MPI
    MPI_Init(&argc, &argv);
    #endif
    
    HDF5Handler h5;

    int domainSize = 64;
    int numIter = 3;
    int testNum = 0;
    double physDomainSize = 1.0;
    std::array<bool, DIM> periodicity;
    periodicity.fill(true);

    InputArgs args;
    
    args.add("domainSize", domainSize);
    args.add("numIter",    numIter);
    args.add("testNum",    testNum);
    args.add("periodic_x", periodicity[0]);
    args.add("periodic_y", periodicity[1]);
#if DIM > 2
    args.add("periodic_z", periodicity[2]);
#endif
    args.parse(argc, argv);
    args.print();

    if (testNum == 0 && procID() == 0)
    {
        std::cout << "TEST 0: MEMTYPE == HOST" << std::endl;
    } else if (testNum == 1 && procID() == 0)
    {
        std::cout << "TEST 1: MEMTYPE == DEVICE" << std::endl;
    } else if (procID() == 0)
    {
        std::cout << "UNKNOWN TEST NUMBER: " << testNum << std::endl;
        return;
    }

    double err[numIter];
    for (int nn = 0; nn < numIter; nn++)
    {
        double dx = physDomainSize / domainSize;
        Box domainBox = Box::Cube(domainSize);
        if (testNum == 0)
        {
            BoxOp_Laplace<double, HOST> op(dx);
            BoxData<double, NUMCOMPS, HOST> Phi(domainBox.grow(op.ghost()));
            BoxData<double, NUMCOMPS, HOST> LPhi(domainBox);
            BoxData<double, NUMCOMPS, HOST> LPhiSln(domainBox);
            BoxData<double, NUMCOMPS, HOST> LPhiErr(domainBox);
            std::array<BoxData<double, NUMCOMPS, HOST>, DIM> Flux;

            forallInPlace_p(sinProd_avg, Phi, dx);
            forallInPlace_p(LsinProd_avg, LPhiSln, dx);
            LPhi.setVal(0);
            LPhiErr.setVal(0);
            
            for (int dir = 0; dir < DIM; dir++)
            {
                Flux[dir].define(domainBox.grow(dir, Side::Hi, 1));
            }

            op(LPhi, Flux, Phi);
            
            LPhi.copyTo(LPhiErr);
            LPhiErr -= LPhiSln;
            
            err[nn] = LPhiErr.absMax(); 
        } else if (testNum == 1)
        {
            BoxOp_Laplace<double, DEVICE> op(dx);
            BoxData<double, NUMCOMPS, DEVICE> Phi(domainBox.grow(op.ghost()));
            BoxData<double, NUMCOMPS, DEVICE> LPhi(domainBox);
            BoxData<double, NUMCOMPS, DEVICE> LPhiSln(domainBox);
            BoxData<double, NUMCOMPS, DEVICE> LPhiErr(domainBox);
            std::array<BoxData<double, NUMCOMPS, DEVICE>, DIM> Flux;

            forallInPlace_p(sinProd_avg, Phi, dx);
            forallInPlace_p(LsinProd_avg, LPhiSln, dx);
            LPhi.setVal(0);
            LPhiErr.setVal(0);
            
            for (int dir = 0; dir < DIM; dir++)
            {
                Flux[dir].define(domainBox.grow(dir, Side::Hi, 1));
            }

            op(LPhi, Flux, Phi);
            
            LPhi.copyTo(LPhiErr);
            LPhiErr -= LPhiSln;
            
            err[nn] = LPhiErr.absMax(); 
        }
        std::cout << "Error: " << err[nn] << std::endl;
        domainSize *= 2;
    }
        
    for (int ii = 1; ii < numIter; ii++)
    {
        std::cout << "Convergence Rate: " << log(err[ii-1] / err[ii]) / log(2.0) << std::endl;
    }

    #ifdef PR_MPI
    MPI_Finalize();
    #endif
    return 0;
}

