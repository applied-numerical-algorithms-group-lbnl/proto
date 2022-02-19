#include "Proto.H"
#include "InterpStencilTP.H"
#include "InputParser.H"
#include "func.H"

using namespace Proto;

int main(int argc, char** argv)
{
    #ifdef PR_MPI
    MPI_Init(&argc, &argv);
    #endif

    // SETUP
    using Proto::pout;
    HDF5Handler h5;
    InputArgs args;
    args.parse();
    args.print();

    int domainSize = 64;
    int boxSize = 16;
    int numIter = 3;
    int ghostSize = 1;
    int refRatio = 2;
    double k = 1;
    double physDomainSize = 1;
    std::array<bool, DIM> periodicity;
    for (int dir = 0; dir < DIM; dir++) { periodicity[dir] = true; }

    args.set("domainSize", &domainSize);
    args.set("boxSize",    &boxSize);
    args.set("numIter",    &numIter);
    args.set("refRatio",   &refRatio);
    args.set("periodic_x", &periodicity[0]);
    args.set("periodic_y", &periodicity[1]);
    
    double err[numIter];
    for (int nn = 0; nn < numIter; nn++)
    {
        double dx = 1.0/domainSize;
        double fdx = dx/refRatio;

        InterpStencilTP<double> interp(4, refRatio);
        std::cout << "span: " << interp.ghost() << std::endl;
        
        Box domainBox = Box::Cube(domainSize).shift(Point::Ones(domainSize));
        Box rangeBox  = domainBox.refine(refRatio);
        BoxData<double> input(domainBox.grow(interp.ghost()));
        BoxData<double> output(rangeBox);
        BoxData<double> solution(rangeBox);
        BoxData<double> error(rangeBox);

        forallInPlace_p(f_wave_avg, input, dx,  k);
        forallInPlace_p(f_wave_avg, solution, fdx, k);

        output.setVal(0);
        
        h5.writePatch(dx, input, "INPUT_N%i", nn);
        interp.apply(output, input);
        h5.writePatch(fdx, output,   "OUTPUT_N%i", nn);
        h5.writePatch(fdx, solution, "SOLUTION_N%i", nn);
        
        output.copyTo(error);
        error -= solution;
        h5.writePatch(fdx, error, "ERROR_N%i", nn);
        err[nn] = error.absMax();

        pout() << "Error: " << err[nn] << std::endl;
        domainSize *= 2;
    }
        
    for (int ii = 1; ii < numIter; ii++)
    {
        pout() << "Convergence Rate: " << log(err[ii-1] / err[ii]) / log(2.0) << std::endl;
    }

    #ifdef PR_MPI
    MPI_Finalize();
    #endif
    return 0;
}

